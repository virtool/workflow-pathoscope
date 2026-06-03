import asyncio
import os
import re
import shlex
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from virtool.caches.utils import derive_key
from virtool.workflow import hooks, step
from virtool.workflow.data.analyses import WFAnalysis
from virtool.workflow.data.cache import CacheHit, WorkflowCache
from virtool.workflow.data.indexes import WFIndex
from virtool.workflow.data.samples import WFSample
from virtool.workflow.data.subtractions import WFSubtraction
from virtool.workflow.runtime.run_subprocess import RunSubprocess
from workflow_pathoscope.utils import (
    run_pathoscope,
    write_default_isolate_fasta,
    write_isolate_fasta,
    write_report,
)

from workflow_pathoscope.rust import (
    find_candidate_otus_with_bowtie2,
    init_logging,
    run_eliminate_subtraction,
)

# Initialize Rust logging to forward to Python logging system
# This enables unified logging across Python and Rust components
init_logging("info")


WORKFLOW_NAME = "pathoscope"
BOWTIE2_BUILD_TOOL = "bowtie2-build"
BOWTIE2_INDEX_SUFFIXES = (
    "1.bt2",
    "2.bt2",
    "3.bt2",
    "4.bt2",
    "rev.1.bt2",
    "rev.2.bt2",
)


async def get_bowtie2_build_version(run_subprocess: RunSubprocess) -> str:
    output = []

    async def collect_stdout(line: bytes) -> None:
        output.append(line.decode())

    await run_subprocess(
        [BOWTIE2_BUILD_TOOL, "--version"],
        stdout_handler=collect_stdout,
    )

    match = re.search(r"\bversion\s+([^\s]+)", "".join(output))

    if match is None:
        raise ValueError("Could not parse bowtie2-build version")

    return match.group(1)


async def get_mapping_index_cache_params(
    index_kind: str,
    parent_id: str,
    run_subprocess: RunSubprocess,
    extra_params: dict[str, str] | None = None,
) -> dict[str, str]:
    tool_version = await get_bowtie2_build_version(run_subprocess)

    params = {
        "index_kind": index_kind,
        "workflow": WORKFLOW_NAME,
        "parent_id": parent_id,
        "tool_name": BOWTIE2_BUILD_TOOL,
        "tool_version": tool_version,
    }

    if extra_params is not None:
        params.update(extra_params)

    return params


async def build_bowtie2_index(
    fasta_path: Path,
    index_prefix: Path,
    proc: int,
    run_subprocess: RunSubprocess,
) -> None:
    await run_subprocess(
        [
            BOWTIE2_BUILD_TOOL,
            "--threads",
            str(proc),
            str(fasta_path),
            str(index_prefix),
        ],
    )


@hooks.on_failure
async def delete_analysis_document(analysis: WFAnalysis):
    await analysis.delete()


@step
async def create_reference_index(
    cache: WorkflowCache,
    index: WFIndex,
    logger,
    proc: int,
    run_subprocess: RunSubprocess,
    reference_index_path: Path,
):
    """Ensure the reference Bowtie2 index exists locally."""
    reference_fasta_path = reference_index_path.parent / "reference.fa"
    reference_index_dir = reference_index_path.parent
    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        run_subprocess,
        {
            "source": "index.json_path",
            "selection": "default_isolates",
        },
    )
    key = derive_key(params)
    log = logger.bind(
        index_kind="reference_mapping_index",
        key=key,
        parent_id=index.id,
        workflow=WORKFLOW_NAME,
    )

    log.info("checking workflow cache")

    result = await cache.get(key, reference_index_dir)

    if isinstance(result, CacheHit):
        log.info("restored cached mapping index", outcome="hit")
        return

    log.info("building mapping index", outcome="miss")

    await asyncio.to_thread(reference_index_dir.mkdir, parents=True, exist_ok=True)

    lengths = await asyncio.to_thread(
        write_default_isolate_fasta,
        index.json_path,
        reference_fasta_path,
    )

    logger.info(
        "assembled default reference fasta",
        count=len(lengths),
        source=str(index.json_path),
    )

    await build_bowtie2_index(
        reference_fasta_path,
        reference_index_path,
        proc,
        run_subprocess,
    )

    await cache.put(key, reference_index_dir, params=params)

    log.info("cached built mapping index", outcome="put")


async def create_subtraction_index(
    cache: WorkflowCache,
    logger,
    proc: int,
    run_subprocess: RunSubprocess,
    subtraction: WFSubtraction,
    work_path: Path,
) -> Path:
    subtraction_index_path = (
        work_path / "subtraction_indexes" / subtraction.id / "subtraction"
    )
    subtraction_index_dir = subtraction_index_path.parent
    params = await get_mapping_index_cache_params(
        "subtraction_mapping_index",
        subtraction.id,
        run_subprocess,
    )
    key = derive_key(params)
    log = logger.bind(
        index_kind="subtraction_mapping_index",
        key=key,
        parent_id=subtraction.id,
        workflow=WORKFLOW_NAME,
    )

    log.info("checking workflow cache")

    result = await cache.get(key, subtraction_index_dir)

    if isinstance(result, CacheHit):
        log.info("restored cached mapping index", outcome="hit")
        return subtraction_index_path

    log.info("building mapping index", outcome="miss")

    await asyncio.to_thread(subtraction_index_dir.mkdir, parents=True, exist_ok=True)

    await build_bowtie2_index(
        subtraction.fasta_path,
        subtraction_index_path,
        proc,
        run_subprocess,
    )

    await cache.put(key, subtraction_index_dir, params=params)

    log.info("cached built mapping index", outcome="put")

    return subtraction_index_path


@step
async def map_default_isolates(
    intermediate: SimpleNamespace,
    logger,
    index: WFIndex,
    proc: int,
    p_score_cutoff: float,
    reference_index_path: Path,
    sample: WFSample,
):
    """Map sample reads to all default isolates to identify candidate OTUs.

    This will be used to identify candidate OTUs.
    """
    logger.info("running bowtie2 directly from rust with streaming")

    candidate_otus = await asyncio.to_thread(
        find_candidate_otus_with_bowtie2,
        str(reference_index_path),
        [str(path) for path in sample.read_paths],
        proc,
        p_score_cutoff,
    )

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

    intermediate.lengths = await asyncio.to_thread(
        write_isolate_fasta,
        {index.get_otu_id_by_sequence_id(id_) for id_ in intermediate.to_otus},
        index.json_path,
        isolate_fasta_path,
    )

    await build_bowtie2_index(
        isolate_fasta_path,
        isolate_index_path,
        proc,
        run_subprocess,
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
    cache: WorkflowCache,
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

    :param cache: workflow cache client
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

    if len(subtractions) == 0:
        logger.info("no subtractions to eliminate reads against")
        # Rename BAM file as no subtraction is needed (saves disk space)
        await asyncio.to_thread(os.rename, isolate_bam_path, subtracted_bam_path)
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
            "processing subtraction",
            id=subtraction.id,
            name=subtraction.name,
        )

        subtraction_index_path = await create_subtraction_index(
            cache,
            logger,
            proc,
            run_subprocess,
            subtraction,
            work_path,
        )

        bowtie_cmd = (
            f"bowtie2 --local --no-unal -N 0 -p {proc} "
            f"-x {shlex.quote(str(subtraction_index_path))} "
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

        eliminated_count = await asyncio.to_thread(
            run_eliminate_subtraction,
            str(current_bam_input_path),
            str(to_subtraction_bam_path),
            str(subtracted_bam_path),
            str(current_fastq_path),
            str(current_fastq_path),
            proc - 1,
        )

        await asyncio.to_thread(to_subtraction_bam_path.unlink)

        current_bam_input_path = await asyncio.to_thread(
            subtracted_bam_path.rename, work_path / "working_isolate.bam"
        )

        subtracted_count += eliminated_count

        logger.info(
            "subtraction complete",
            id=subtraction.id,
            name=subtraction.name,
            eliminated_this_subtraction=eliminated_count,
            total_eliminated=subtracted_count,
        )

    await asyncio.to_thread(current_bam_input_path.rename, subtracted_bam_path)

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
        run_pathoscope, subtracted_bam_path, p_score_cutoff
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

    hits = []

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
