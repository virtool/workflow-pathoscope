import asyncio
import shlex
import shutil
from collections import defaultdict
from logging import getLogger
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, TextIO, Set

import rust_utils
import aiofiles
import aiofiles.os
from virtool_workflow import hooks, step
from virtool_workflow.data_model.analysis import WFAnalysis
from virtool_workflow.data_model.indexes import WFIndex
from virtool_workflow.data_model.subtractions import WFSubtraction
from virtool_workflow.runtime.run_subprocess import RunSubprocess

from pathoscope import SamLine, calculate_coverage
from pathoscope import run as run_pathoscope
from pathoscope import write_report

logger = getLogger("pathoscope")

BAD_FIRST_SAM_CHARACTERS = {"\n", "@", "#"}


def read_fastq_grouped_lines(fastq_file: TextIO):
    while True:
        fastq_read = (
            fastq_file.readline(),
            fastq_file.readline(),
            fastq_file.readline(),
            fastq_file.readline(),
        )

        if "" in fastq_read:
            return

        yield fastq_read


def subtract_fastq(
    current_fastq_path: Path, new_fastq_path: Path, subtracted_reads: Set[str]
):
    with open(current_fastq_path, "r") as current_fastq_file, open(
        new_fastq_path, "w"
    ) as new_fastq_file:
        for record in read_fastq_grouped_lines(current_fastq_file):
            if record[0].strip("@\n") not in subtracted_reads:
                new_fastq_file.write("".join(record))


@hooks.on_failure
async def delete_analysis_document(analysis_provider):
    await analysis_provider.delete()


@hooks.on_result
async def upload_results(results, analysis_provider):
    await analysis_provider.upload_result(results)


@step
async def map_default_isolates(
    intermediate: SimpleNamespace,
    read_file_names: str,
    index: WFIndex,
    proc: int,
    p_score_cutoff: float,
    run_subprocess: RunSubprocess,
):
    """
    Map sample reads to all default isolates to identify candidate OTUs.

    This will be used to identify candidate OTUs.
    """

    async def stdout_handler(line: bytes):
        line = line.decode()

        if line[0] == "#" or line[0] == "@":
            return

        sam_line = SamLine(line)

        if sam_line.unmapped:
            return

        if sam_line.ref_id == "*":
            return

        if sam_line.score < p_score_cutoff:
            return

        intermediate.to_otus.add(sam_line.ref_id)

    await run_subprocess(
        [
            "bowtie2",
            "-p",
            str(proc),
            "--no-unal",
            "--local",
            "--score-min",
            "L,20,1.0",
            "-N",
            "0",
            "-L",
            "15",
            "-x",
            str(index.bowtie_path),
            "-U",
            read_file_names,
        ],
        stdout_handler=stdout_handler,
    )

    logger.info("Found %i potential OTUs.", len(intermediate.to_otus))


@step
async def build_isolate_index(
    index: WFIndex,
    intermediate: SimpleNamespace,
    isolate_fasta_path: Path,
    isolate_index_path: Path,
    run_subprocess: RunSubprocess,
    proc: int,
):
    """
    Build a mapping index containing all isolates of candidate OTUs.
    """
    intermediate.lengths = await index.write_isolate_fasta(
        [index.get_otu_id_by_sequence_id(id_) for id_ in intermediate.to_otus],
        isolate_fasta_path,
        proc,
    )

    await run_subprocess(
        [
            "bowtie2-build",
            "--threads",
            str(proc),
            str(isolate_fasta_path),
            str(isolate_index_path),
        ]
    )


@step
async def map_isolates(
    read_file_names: str,
    intermediate: SimpleNamespace,
    isolate_fastq_path: Path,
    isolate_index_path: Path,
    isolate_sam_path: Path,
    run_subprocess: RunSubprocess,
    proc: int,
    p_score_cutoff: float,
):
    """Map sample reads to the all isolate index."""
    isolate_high_scores = defaultdict(float)

    async with aiofiles.open(isolate_sam_path, "w") as f:

        async def stdout_handler(line: bytes):
            line = line.decode()

            if not line or line[0] in BAD_FIRST_SAM_CHARACTERS:
                return

            sam_line = SamLine(line)

            if sam_line.unmapped:
                return

            if sam_line.ref_id == "*":
                return

            # Skip if the p_score does not meet the minimum cutoff.
            if sam_line.score < p_score_cutoff:
                return

            if sam_line.score > isolate_high_scores[sam_line.read_id]:
                isolate_high_scores[sam_line.read_id] = sam_line.score

            await f.write(f"{line}")

        command = [
            "bowtie2",
            "-p",
            str(proc),
            "--no-hd",
            "--no-sq",
            "--no-unal",
            "--local",
            "--score-min",
            "L,20,1.0",
            "-N",
            "0",
            "-L",
            "15",
            "-k",
            "100",
            "--al",
            str(isolate_fastq_path),
            "-x",
            str(isolate_index_path),
            "-U",
            read_file_names,
        ]

        await run_subprocess(command, stdout_handler=stdout_handler)

    intermediate.isolate_high_scores = dict(isolate_high_scores)


@step
async def eliminate_subtraction(
    isolate_fastq_path: Path,
    isolate_sam_path: Path,
    proc: int,
    results: Dict[str, Any],
    run_subprocess: RunSubprocess,
    subtractions: List[WFSubtraction],
    subtracted_sam_path: Path,
    work_path: Path,
):
    """
    Remove reads that map better to a subtraction than to a reference.

    The input to this step is the reads that aligned to an isolate at least once in the
    previous step. We will align these against a subtraction (plant) genome. If the
    alignment score is higher against the subtraction, we drop alignments involving the
    read from the SAM from the previous step and write the reduced one to
    `subtracted_sam_path`.

    TODO: Go through all `subtractions` instead of just the first `subtraction`.
          Efficiently eliminate from the subtraction SAM any reads that map to any of
          subtractions.

    :param isolate_fastq_path: path to the FASTQ file containing reads that aligned to the isolates
    :param isolate_sam_path: path to the SAM file of alignments to the isolates
    :param proc: number of processors to use
    :param results: the results to send to the api when the workflow is complete
    :param run_subprocess: runs a subprocess with error handling
    :param subtractions: the subtraction to align and eliminate reads against
    :param subtracted_sam_path: path to the SAM file with subtraction-mapped reads removed
    :param work_path: path to the workflow working directory
    """

    if len(subtractions) == 0:
        logger.info("No subtractions to eliminate reads against. Skipping step.")
        await asyncio.to_thread(shutil.copyfile, isolate_sam_path, subtracted_sam_path)
        results["subtracted_count"] = 0
        return

    current_fastq_path = work_path / "current_fastq.fq"
    to_subtraction_sam_path = work_path / "to_subtraction.sam"

    # copy the original fastq file into a working fastq file
    # as to not disrupt possible uses elsewhere
    await asyncio.to_thread(shutil.copyfile, isolate_fastq_path, current_fastq_path)

    # The SAM file that reads should be subtracted from if they map better to a
    # subtraction.
    current_sam_input_path = isolate_sam_path

    subtracted_count = 0

    for subtraction in subtractions:
        # Map reads to the subtraction.
        await run_subprocess(
            [
                "bowtie2",
                "--local",
                "--no-unal",
                "--no-hd",
                "--no-sq",
                "-N",
                "0",
                "-p",
                str(proc),
                "-x",
                shlex.quote(str(subtraction.bowtie2_index_path)),
                "-U",
                str(current_fastq_path),
                "-S",
                str(to_subtraction_sam_path),
            ]
        )

        rust_utils.run_eliminate_subtraction(
            str(current_sam_input_path),
            str(to_subtraction_sam_path),
            str(subtracted_sam_path),
        )

        await aiofiles.os.remove(to_subtraction_sam_path)

        current_sam_input_path = work_path / "working_isolate.sam"

        await asyncio.to_thread(
            shutil.copyfile, subtracted_sam_path, current_sam_input_path
        )

        async with aiofiles.open("subtracted_read_ids.txt", "r") as f:
            subtracted_reads = {str(line).strip("@\n") async for line in f}

        subtracted_count += len(subtracted_reads)

        # Rewrite the input FASTQ file to exclude reads that mapped better to a
        # subtraction.
        new_fastq_path = work_path / "new_fastq.fq"

        await asyncio.to_thread(
            subtract_fastq, current_fastq_path, new_fastq_path, subtracted_reads
        )

        # Overwrite the previous input FASTA with the subtracted one.
        await asyncio.to_thread(shutil.copyfile, new_fastq_path, current_fastq_path)

        logger.info(
            "Subtracted reads that mapped better to a subtraction subtraction_id=%s subtraction_name=%s count=%i",
            subtraction.id,
            subtraction.name,
            subtracted_count,
        )

    results["subtracted_count"] = subtracted_count


@step
async def reassignment(
    analysis: WFAnalysis,
    index: WFIndex,
    intermediate: SimpleNamespace,
    p_score_cutoff: float,
    results,
    subtracted_sam_path: Path,
    work_path: Path,
):
    """
    Run the Pathoscope reassignment algorithm.

    Tab-separated output is written to ``pathoscope.tsv``. The results are also parsed
    and saved to `intermediate.coverage`.
    """
    reassigned_path = work_path / "reassigned.sam"

    logger.info(
        "Running Pathoscope subtracted_sam_path=%s reassigned_path=%s",
        subtracted_sam_path,
        reassigned_path,
    )

    (
        best_hit_initial_reads,
        best_hit_initial,
        level_1_initial,
        level_2_initial,
        best_hit_final_reads,
        best_hit_final,
        level_1_final,
        level_2_final,
        init_pi,
        pi,
        refs,
        reads,
    ) = await asyncio.to_thread(
        run_pathoscope, subtracted_sam_path, reassigned_path, p_score_cutoff
    )

    read_count = len(reads)

    report_path = work_path / "report.tsv"

    logger.info("Writing report report_path=%s", report_path)

    report = await asyncio.to_thread(
        write_report,
        report_path,
        pi,
        refs,
        read_count,
        init_pi,
        best_hit_initial,
        best_hit_initial_reads,
        best_hit_final,
        best_hit_final_reads,
        level_1_initial,
        level_2_initial,
        level_1_final,
        level_2_final,
    )

    analysis.upload(report_path, "tsv")

    logger.info("Calculating coverage")

    intermediate.coverage = await asyncio.to_thread(
        calculate_coverage, reassigned_path, intermediate.lengths
    )

    logger.info("Preparing hits")

    hits = list()

    for sequence_id, hit in report.items():
        otu_id = index.get_otu_id_by_sequence_id(sequence_id)

        hit["id"] = sequence_id

        # Attach "otu" (id, version) to the hit.
        hit["otu"] = {"id": otu_id, "version": index.manifest[otu_id]}

        # Get the coverage for the sequence.
        hit_coverage = intermediate.coverage[sequence_id]

        hit["align"] = hit_coverage

        # Calculate coverage and attach to hit.
        hit["coverage"] = round(1 - hit_coverage.count(0) / len(hit_coverage), 3)

        # Calculate depth and attach to hit.
        hit["depth"] = round(sum(hit_coverage) / len(hit_coverage))

        hits.append(hit)

    results.update({"read_count": read_count, "hits": hits})

    intermediate.reassigned_path = reassigned_path
    intermediate.report_path = report_path
