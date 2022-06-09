import shlex
from collections import defaultdict
from logging import getLogger
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import aiofiles
import aiofiles.os
from virtool_workflow import fixture, hooks, step
from virtool_workflow.analysis.analysis import Analysis
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.subtractions import Subtraction
from virtool_workflow.execution.run_subprocess import RunSubprocess

from pathoscope import SamLine, calculate_coverage
from pathoscope import run as run_pathoscope
from pathoscope import write_report

logger = getLogger(__name__)


@fixture
def index(indexes: List[Index]):
    return indexes[0]


@fixture
def subtraction(subtractions: List[Subtraction]):
    return subtractions[0]


@fixture
def intermediate():
    """A namespace for storing intermediate values."""
    return SimpleNamespace(
        isolate_high_scores={},
        to_otus=set(),
    )


@fixture
def isolate_path(work_path: Path):
    path = work_path / "isolates"
    path.mkdir()

    return path


@fixture
def isolate_fasta_path(isolate_path: Path):
    return isolate_path / "isolate_index.fa"


@fixture
def isolate_fastq_path(isolate_path: Path):
    return isolate_path / "isolate_mapped.fq"


@fixture
def isolate_index_path(isolate_path: Path):
    return isolate_path / "isolates"


@fixture
def isolate_sam_path(isolate_path: Path):
    return isolate_path / "to_isolates.sam"


@fixture
def p_score_cutoff():
    return 0.01


@fixture
def read_file_names(sample) -> str:
    return ",".join(str(path) for path in sample.read_paths)


@fixture
def reassigned_sam_path(work_path: Path):
    """The path to the SAM file after Pathoscope reassignment."""
    return work_path / "reassigned.sam"


@fixture
def subtracted_sam_path(work_path: Path) -> Path:
    """The path to the SAM file after subtraction reads have been eliminated."""
    return work_path / "subtracted.sam"


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
    index: Index,
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

    logger.info(f"Found {len(intermediate.to_otus)} potential OTUs.")


@step
async def build_isolate_index(
    index: Index,
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
        ],
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

            if line[0] == "@" or line == "#":
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

            await f.write(f"{line}\n")

        command = [
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
    subtraction: Subtraction,
    subtracted_sam_path: Path,
    work_path: Path,
):
    """
    Remove reads that map better to a subtraction than to a reference.
    """
    to_subtraction_sam_path = work_path / "to_subtraction.sam"

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
            str(isolate_fastq_path),
            "-S",
            str(to_subtraction_sam_path),
        ]
    )

    await run_subprocess(
        [
            "./eliminate_subtraction",
            str(isolate_sam_path),
            str(to_subtraction_sam_path),
            str(subtracted_sam_path),
        ]
    )

    await aiofiles.os.remove(to_subtraction_sam_path)

    subtracted_count = 0

    async with aiofiles.open("subtracted_read_ids.txt", "r") as f:
        async for line in f:
            if line:
                subtracted_count += 1

    logger.info(
        f"Subtracted {subtracted_count} reads that mapped better to a subtraction."
    )

    results["subtracted_count"] = subtracted_count


@step
async def reassignment(
    analysis: Analysis,
    index: Index,
    intermediate: SimpleNamespace,
    p_score_cutoff: float,
    results,
    run_in_executor,
    subtracted_sam_path: Path,
    work_path: Path,
):
    """
    Run the Pathoscope reassignment algorithm.

    Tab-separated output is written to ``pathoscope.tsv``. The results are also parsed
    and saved to `intermediate.coverage`.

    """
    reassigned_path = work_path / "reassigned.sam"

    logger.info(f"subtracted_sam_path: {subtracted_sam_path}")
    logger.info(f"reassigned_path: {reassigned_path}")

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
    ) = await run_in_executor(
        run_pathoscope,
        subtracted_sam_path,
        reassigned_path,
        p_score_cutoff,
    )

    read_count = len(reads)

    report_path = work_path / "report.tsv"

    report = await run_in_executor(
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

    intermediate.coverage = await run_in_executor(
        calculate_coverage, reassigned_path, intermediate.lengths
    )

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
