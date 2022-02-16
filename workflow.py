import shlex
from pathlib import Path
from types import SimpleNamespace
from typing import List

import aiofiles
from virtool_workflow import fixture, hooks, step
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.subtractions import Subtraction
from virtool_workflow.execution.run_subprocess import RunSubprocess

from pathoscope import calculate_coverage, find_sam_align_score, \
    replace_after_subtraction, run as run_pathoscope, subtract, write_report


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
def isolate_vta_path(isolate_path: Path):
    return isolate_path / "to_isolates.vta"


@fixture
def subtraction_vta_path(work_path: Path):
    return work_path / "subtracted.vta"


@fixture
def p_score_cutoff():
    return 0.01


@fixture
def read_file_names(sample) -> str:
    return ",".join(str(path) for path in sample.read_paths)


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
        run_subprocess: RunSubprocess
):
    """
    Map sample reads to all default isolates to identify candidate OTUs.

    This will be used to identify candidate OTUs.

    """
    async def stdout_handler(line: bytes):
        line = line.decode()

        if line[0] == "#" or line[0] == "@":
            return

        fields = line.split("\t")

        # Bitwise FLAG - 0x4: segment unmapped
        if int(fields[1]) & 0x4 == 4:
            return

        ref_id = fields[2]

        if ref_id == "*":
            return

        # Skip if the p_score does not meet the minimum cutoff.
        if find_sam_align_score(fields) < p_score_cutoff:
            return

        intermediate.to_otus.add(ref_id)

    await run_subprocess(
        [
            "bowtie2",
            "-p", str(proc),
            "--no-unal",
            "--local",
            "--score-min", "L,20,1.0",
            "-N", "0",
            "-L", "15",
            "-x", str(index.bowtie_path),
            "-U", read_file_names,
        ],
        stdout_handler=stdout_handler
    )


@step
async def build_isolate_index(
        index: Index,
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
        proc
    )

    await run_subprocess(
        [
            "bowtie2-build",
            "--threads", str(proc),
            str(isolate_fasta_path),
            str(isolate_index_path)
        ],
    )

    return "Built isolate index."


@step
async def map_isolates(
        read_file_names: str,
        isolate_fastq_path: Path,
        isolate_index_path: Path,
        isolate_vta_path: Path,
        run_subprocess: RunSubprocess,
        proc: int,
        p_score_cutoff: float
):
    """Map sample reads to the all isolate index."""
    async with aiofiles.open(isolate_vta_path, "w") as f:
        async def stdout_handler(line: bytes):
            line = line.decode()

            if line[0] == "@" or line == "#":
                return

            fields = line.split("\t")

            # Bitwise FLAG - 0x4 : segment unmapped
            if int(fields[1]) & 0x4 == 4:
                return

            ref_id = fields[2]

            if ref_id == "*":
                return

            p_score = find_sam_align_score(fields)

            # Skip if the p_score does not meet the minimum cutoff.
            if p_score < p_score_cutoff:
                return

            await f.write(",".join([
                fields[0],  # read_id
                ref_id,
                fields[3],  # pos
                str(len(fields[9])),  # length
                str(p_score)
            ]) + "\n")

        command = [
            "bowtie2",
            "-p", str(proc),
            "--no-unal",
            "--local",
            "--score-min", "L,20,1.0",
            "-N", "0",
            "-L", "15",
            "-k", "100",
            "--al", str(isolate_fastq_path),
            "-x", str(isolate_index_path),
            "-U", read_file_names
        ]

        await run_subprocess(command, stdout_handler=stdout_handler)

    return "Mapped sample reads to isolate index."


@step
async def map_subtractions(
        intermediate: SimpleNamespace,
        subtraction: Subtraction,
        isolate_fastq_path,
        run_subprocess: RunSubprocess,
        proc: int,
):
    """
    Map isolate-mapped reads to the subtraction.
    """
    to_subtraction = {}

    async def stdout_handler(line: bytes):
        line = line.decode()

        if line[0] == "@" or line == "#":
            return

        fields = line.split("\t")

        # Bitwise FLAG - 0x4 : segment unmapped
        if int(fields[1]) & 0x4 == 4:
            return

        # No ref_id assigned.
        if fields[2] == "*":
            return

        to_subtraction[fields[0]] = find_sam_align_score(fields)

    command = [
        "bowtie2",
        "--local",
        "-N", "0",
        "-p", str(proc),
        "-x", shlex.quote(str(subtraction.bowtie2_index_path)),
        "-U", str(isolate_fastq_path)
    ]

    await run_subprocess(command, stdout_handler=stdout_handler)

    intermediate.to_subtraction = to_subtraction


@step
async def subtract_mapping(
        intermediate: SimpleNamespace,
        isolate_vta_path: Path,
        results: dict,
        run_in_executor,
        work_path: Path
):
    """Discard reads that map better to the subtraction than an OTU."""
    output_path = work_path / "subtracted.vta"

    subtracted_count = await subtract(
        isolate_vta_path,
        output_path,
        intermediate.to_subtraction
    )

    await run_in_executor(
        replace_after_subtraction,
        output_path,
        isolate_vta_path
    )

    results["subtracted_count"] = subtracted_count


@step
async def reassignment(
        intermediate,
        results,
        run_in_executor,
        isolate_vta_path: Path,
        index: Index,
        p_score_cutoff: float,
        work_path: Path
):
    """
    Run the Pathoscope reassignment algorithm.

    Tab-separated output is written to ``pathoscope.tsv``. The results are also parsed
    and saved to `intermediate.coverage`.
    """
    reassigned_path = work_path / "reassigned.vta"

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
        reads
    ) = await run_in_executor(
        run_pathoscope,
        isolate_vta_path,
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
        level_2_final
    )

    intermediate.coverage = await run_in_executor(
        calculate_coverage,
        reassigned_path,
        intermediate.lengths
    )

    hits = list()

    for sequence_id, hit in report.items():
        otu_id = index.get_otu_id_by_sequence_id(sequence_id)

        hit["id"] = sequence_id

        # Attach "otu" (id, version) to the hit.
        hit["otu"] = {
            "id": otu_id,
            "version": index.manifest[otu_id]
        }

        # Get the coverage for the sequence.
        hit_coverage = intermediate.coverage[sequence_id]

        hit["align"] = hit_coverage

        # Calculate coverage and attach to hit.
        hit["coverage"] = round(
            1 - hit_coverage.count(0) / len(hit_coverage), 3)

        # Calculate depth and attach to hit.
        hit["depth"] = round(sum(hit_coverage) / len(hit_coverage))

        hits.append(hit)

    results.update({
        "read_count": read_count,
        "hits": hits
    })

    intermediate.reassigned_path = reassigned_path
    intermediate.report_path = report_path
