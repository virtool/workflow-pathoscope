import asyncio
import shlex
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import aiofiles
import aiofiles.os
from virtool_workflow import hooks, step
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.indexes import WFIndex
from virtool_workflow.data.samples import WFSample
from virtool_workflow.data.subtractions import WFSubtraction
from virtool_workflow.runtime.run_subprocess import RunSubprocess

from workflow_pathoscope.utils import (
    run_pathoscope,
    write_report,
)
from workflow_pathoscope.rust import (
    parse_isolate_scores,
    find_candidate_otus_from_bytes,
    eliminate_subtraction_and_filter_fastq,
    init_logging,
)

# Initialize Rust logging to forward to Python logging system
# This enables unified logging across Python and Rust components
init_logging("info")


@hooks.on_failure
async def delete_analysis_document(analysis: WFAnalysis):
    await analysis.delete()


@step
async def map_default_isolates(
    intermediate: SimpleNamespace,
    logger,
    index: WFIndex,
    proc: int,
    p_score_cutoff: float,
    sample: WFSample,
):
    """Map sample reads to all default isolates to identify candidate OTUs.

    This will be used to identify candidate OTUs.
    """
    # Run bowtie2 and capture output in memory
    command = [
        "bowtie2",
        "-p",
        str(proc),
        "--local",
        "--no-unal",
        "--score-min",
        "L,20,1.0",
        "-N",
        "0",
        "-L",
        "15",
        "-x",
        index.bowtie_path,
        "-U",
        ",".join(str(path) for path in sample.read_paths),
    ]

    logger.info("running bowtie2 and capturing output")

    # Execute bowtie2 and capture stdout
    process = await asyncio.create_subprocess_exec(
        *command,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )

    stdout, stderr = await process.communicate()

    if process.returncode != 0:
        logger.error(
            "bowtie2 failed", returncode=process.returncode, stderr=stderr.decode()
        )
        raise RuntimeError(f"bowtie2 failed with return code {process.returncode}")

    # Use Rust implementation to process the SAM bytes directly
    logger.info("extracting candidate otus from memory")

    candidate_otus = await asyncio.to_thread(
        find_candidate_otus_from_bytes,
        stdout,
        p_score_cutoff,
    )

    # Convert Rust HashSet to Python set
    intermediate.to_otus = set(candidate_otus)

    logger.info("found candidate otus", count=len(intermediate.to_otus))


@step
async def build_isolate_index(
    index: WFIndex,
    intermediate: SimpleNamespace,
    isolate_fasta_path: Path,
    isolate_index_path: Path,
    run_subprocess: RunSubprocess,
    proc: int,
):
    """Build a mapping index containing all isolates of candidate OTUs."""
    intermediate.lengths = await index.write_isolate_fasta(
        [index.get_otu_id_by_sequence_id(id_) for id_ in intermediate.to_otus],
        isolate_fasta_path,
    )

    await run_subprocess(
        [
            "bowtie2-build",
            "--threads",
            str(proc),
            str(isolate_fasta_path),
            str(isolate_index_path),
        ],
    )


@step
async def map_isolates(
    isolate_fastq_path: Path,
    isolate_index_path: Path,
    isolate_bam_path: Path,
    proc: int,
    run_subprocess: RunSubprocess,
    sample: WFSample,
):
    """Map sample reads to the all isolate index."""
    read_paths = ",".join(str(path) for path in sample.read_paths)

    bowtie_cmd = (
        f"bowtie2 -p {proc} --no-unal --local --score-min L,20,1.0 -N 0 -L 15 -k 100 "
        f"--al {isolate_fastq_path} -x {isolate_index_path} -U {read_paths}"
    )

    samtools_cmd = f"samtools view -bS - -o {isolate_bam_path}"

    await run_subprocess(
        [
            "bash",
            "-c",
            f"{bowtie_cmd} | {samtools_cmd}",
        ]
    )


@step
async def eliminate_subtraction(
    intermediate: SimpleNamespace,
    isolate_fastq_path: Path,
    isolate_bam_path: Path,
    logger,
    p_score_cutoff: float,
    proc: int,
    results: dict[str, Any],
    run_subprocess: RunSubprocess,
    subtractions: list[WFSubtraction],
    subtracted_bam_path: Path,
    work_path: Path,
):
    """Remove reads that map better to a subtraction than to a reference.

    The input to this step is the reads that aligned to an isolate at least once in the
    previous step. We will align these against a subtraction (plant) genome. If the
    alignment score is higher against the subtraction, we drop alignments involving the
    read from the BAM from the previous step and write the reduced one to
    `subtracted_bam_path`.

    :param intermediate: intermediate data storage for the workflow
    :param isolate_fastq_path: path to the FASTQ file containing reads that aligned to the isolates
    :param isolate_bam_path: path to the BAM file of alignments to the isolates
    :param logger: workflow logger
    :param p_score_cutoff: minimum p_score cutoff for alignments
    :param proc: number of processors to use
    :param results: the results to send to the api when the workflow is complete
    :param run_subprocess: runs a subprocess with error handling
    :param subtractions: the subtraction to align and eliminate reads against
    :param subtracted_bam_path: path to the BAM file with subtraction-mapped reads removed
    :param work_path: path to the workflow working directory
    """
    # Parse isolate scores from the BAM file
    logger.info("parsing isolate scores from bam file")
    intermediate.isolate_high_scores = await asyncio.to_thread(
        parse_isolate_scores,
        str(isolate_bam_path),
        p_score_cutoff,
    )
    logger.info("found isolate scores", count=len(intermediate.isolate_high_scores))

    if len(subtractions) == 0:
        logger.info("No subtractions to eliminate reads against. Skipping step.")
        # Rename BAM file as no subtraction is needed (saves disk space)
        await aiofiles.os.rename(isolate_bam_path, subtracted_bam_path)
        results["subtracted_count"] = 0
        return

    current_fastq_path = work_path / "current_fastq.fq"
    to_subtraction_bam_path = work_path / "to_subtraction.bam"

    # copy the original fastq file into a working fastq file
    # as to not disrupt possible uses elsewhere
    await asyncio.to_thread(shutil.copyfile, isolate_fastq_path, current_fastq_path)

    # The file that reads should be subtracted from if they map better to a
    # subtraction. Start with BAM, then use working BAM for subsequent iterations.
    current_bam_input_path = isolate_bam_path

    subtracted_count = 0

    for subtraction in subtractions:
        logger.info(
            "Processing subtraction",
            id=subtraction.id,
            name=subtraction.name,
        )

        bowtie_cmd = (
            f"bowtie2 --local --no-unal -N 0 -p {proc} "
            f"-x {shlex.quote(str(subtraction.bowtie2_index_path))} "
            f"-U {current_fastq_path}"
        )

        samtools_cmd = f"samtools view -bS - -o {to_subtraction_bam_path}"

        await run_subprocess(
            [
                "bash",
                "-c",
                f"{bowtie_cmd} | {samtools_cmd}",
            ],
        )

        # Combined operation: eliminate subtraction reads from BAM and filter FASTQ
        eliminated_count = await asyncio.to_thread(
            eliminate_subtraction_and_filter_fastq,
            str(current_bam_input_path),
            str(to_subtraction_bam_path),
            str(subtracted_bam_path),
            str(current_fastq_path),
            str(current_fastq_path),  # Write directly back to current_fastq_path
        )

        await aiofiles.os.remove(to_subtraction_bam_path)

        current_bam_input_path = work_path / "working_isolate.bam"

        await aiofiles.os.rename(subtracted_bam_path, current_bam_input_path)

        subtracted_count += eliminated_count

        logger.info(
            "Completed subtraction - reads eliminated",
            id=subtraction.id,
            name=subtraction.name,
            eliminated_this_subtraction=eliminated_count,
            total_eliminated=subtracted_count,
        )

    # Rename final working file back to expected output path
    await aiofiles.os.rename(current_bam_input_path, subtracted_bam_path)

    results["subtracted_count"] = subtracted_count


@step
async def reassignment(
    analysis: WFAnalysis,
    index: WFIndex,
    intermediate: SimpleNamespace,
    logger,
    p_score_cutoff: float,
    results,
    subtracted_bam_path: Path,
    work_path: Path,
):
    """Run the Pathoscope reassignment algorithm.

    Tab-separated output is written to ``pathoscope.tsv``. The results are also parsed
    and saved to `intermediate.coverage`.
    """
    logger.info(
        "running pathoscope",
        subtracted_path=subtracted_bam_path,
    )

    pathoscope_results = await asyncio.to_thread(
        run_pathoscope,
        subtracted_bam_path,
        p_score_cutoff,
        intermediate.lengths,
    )

    report_path = work_path / "report.tsv"

    logger.info("writing report", path=report_path)

    report = await asyncio.to_thread(
        write_report,
        report_path,
        pathoscope_results,
    )

    await analysis.upload_file(report_path, "tsv")

    logger.info("preparing hits")

    hits = list()

    for sequence_id, hit in report.items():
        otu_id = index.get_otu_id_by_sequence_id(sequence_id)

        hit["id"] = sequence_id

        # Attach "otu" (id, version) to the hit.
        hit["otu"] = {"id": otu_id, "version": index.manifest[otu_id]}

        # Get the coverage for the sequence.
        hit_coverage = pathoscope_results.coverage[sequence_id]

        hit["align"] = hit_coverage

        # Calculate coverage and attach to hit.
        hit["coverage"] = round(1 - hit_coverage.count(0) / len(hit_coverage), 3)

        # Calculate depth and attach to hit.
        hit["depth"] = round(sum(hit_coverage) / len(hit_coverage))

        hits.append(hit)

    await analysis.upload_result(
        {**results, "read_count": len(pathoscope_results.reads), "hits": hits}
    )
