import asyncio
import csv
import gzip
import json
import re
from pathlib import Path

from virtool.caches.utils import derive_key
from virtool.workflow.data.cache import CacheHit, WorkflowCache
from virtool.workflow.runtime.run_subprocess import RunSubprocess
from virtool.workflow.utils import get_workflow_version

from workflow_pathoscope.rust import PathoscopeResults, run_expectation_maximization


BOWTIE2_BUILD_TOOL = "bowtie2-build"
BOWTIE2_BUILD_VERSION_ATTR = "_bowtie2_build_version"
WORKFLOW_NAME = "pathoscope"
WORKFLOW_VERSION = get_workflow_version()


async def get_bowtie2_build_version(run_subprocess: RunSubprocess) -> str:
    cached_version = getattr(run_subprocess, BOWTIE2_BUILD_VERSION_ATTR, None)

    if cached_version is not None:
        return cached_version

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

    version = match.group(1)
    setattr(run_subprocess, BOWTIE2_BUILD_VERSION_ATTR, version)

    return version


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
        "workflow_version": WORKFLOW_VERSION,
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


async def create_mapping_index(
    cache: WorkflowCache,
    logger,
    proc: int,
    run_subprocess: RunSubprocess,
    *,
    fasta_path: Path,
    index_kind: str,
    index_prefix: Path,
    parent_id: str,
    extra_params: dict[str, str] | None = None,
) -> None:
    index_dir = index_prefix.parent
    cache_restore_parent = index_dir.parent
    params = await get_mapping_index_cache_params(
        index_kind,
        parent_id,
        run_subprocess,
        extra_params,
    )
    key = derive_key(params)
    log = logger.bind(
        index_kind=index_kind,
        key=key,
        parent_id=parent_id,
        workflow=WORKFLOW_NAME,
    )

    log.info("checking workflow cache")

    result = await cache.get(key, cache_restore_parent)

    if isinstance(result, CacheHit):
        log.info("restored cached mapping index", outcome="hit")
        return

    log.info("building mapping index", outcome="miss")

    await asyncio.to_thread(index_dir.mkdir, parents=True, exist_ok=True)

    await build_bowtie2_index(
        fasta_path,
        index_prefix,
        proc,
        run_subprocess,
    )

    created = await cache.put(key, index_dir, params=params)

    if created:
        log.info("cached built mapping index", outcome="put")
    else:
        log.info("mapping index cache already exists", outcome="put_skipped")


def _open_json(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")

    return open(path)


def _get_reference_otus(reference_data):
    if isinstance(reference_data, dict):
        return reference_data["otus"]

    return reference_data


def write_default_isolate_fasta(
    json_path: Path,
    target_path: Path,
) -> dict[str, int]:
    """Generate a FASTA file containing only default isolates from a reference JSON."""
    lengths = {}

    with _open_json(json_path) as f_json, open(target_path, "w") as f_target:
        for otu in _get_reference_otus(json.load(f_json)):
            for isolate in otu["isolates"]:
                if not isolate["default"]:
                    continue

                for sequence in isolate["sequences"]:
                    f_target.write(f">{sequence['_id']}\n{sequence['sequence']}\n")
                    lengths[sequence["_id"]] = len(sequence["sequence"])

    return lengths


def write_isolate_fasta(
    otu_ids: set[str],
    json_path: Path,
    target_path: Path,
) -> dict[str, int]:
    """Generate a FASTA file for all the isolates of the OTUs specified by ``otu_ids``.

    :param otu_ids: the list of OTU IDs for which to generate and index
    :param json_path: the path to the reference index json file
    :param target_path: the path to write the fasta file to
    :return: a dictionary of the lengths of all sequences keyed by their IDS

    """
    lengths = {}

    with _open_json(json_path) as f_json, open(target_path, "w") as f_target:
        for otu in _get_reference_otus(json.load(f_json)):
            if otu["_id"] in otu_ids:
                for isolate in otu["isolates"]:
                    for sequence in isolate["sequences"]:
                        f_target.write(f">{sequence['_id']}\n{sequence['sequence']}\n")
                        lengths[sequence["_id"]] = len(sequence["sequence"])

    return lengths


def write_report(
    path,
    pathoscope_results: PathoscopeResults,
):
    """Write pathoscope results to TSV report and return processed results dict.

    Float values are rounded to 10 decimal places for consistent precision
    across Rust/Python computations.
    """
    read_count = len(pathoscope_results.reads)

    tmp = zip(
        pathoscope_results.pi,
        pathoscope_results.refs,
        pathoscope_results.init_pi,
        pathoscope_results.best_hit_initial,
        pathoscope_results.best_hit_initial_reads,
        pathoscope_results.best_hit_final,
        pathoscope_results.best_hit_final_reads,
        pathoscope_results.level_1_initial,
        pathoscope_results.level_2_initial,
        pathoscope_results.level_1_final,
        pathoscope_results.level_2_final,
    )

    tmp = sorted(tmp, reverse=True)

    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = zip(*tmp)

    end = 0

    for i, _ in enumerate(x10):
        if x1[i] < 0.01 and x10[i] <= 0 and x11[i] <= 0:
            break

        end += 1

    with open(path, "w") as handle:
        csv_writer = csv.writer(handle, delimiter="\t")

        csv_writer.writerow(
            [
                "Total Number of Aligned Reads:",
                read_count,
                "Total Number of Mapped Genomes:",
                len(pathoscope_results.refs),
            ],
        )

        csv_writer.writerow(
            [
                "Genome",
                "Final Guess",
                "Final Best Hit",
                "Final Best Hit Read Numbers",
                "Final High Confidence Hits",
                "Final Low Confidence Hits",
                "Initial Guess",
                "Initial Best Hit",
                "Initial Best Hit Read Numbers",
                "Initial High Confidence Hits",
                "Initial Low Confidence Hits",
            ],
        )

        # Change the column order with zip.
        csv_writer.writerows(
            zip(
                x2[:end],
                x1[:end],
                x6[:end],
                x7[:end],
                x10[:end],
                x11[:end],
                x3[:end],
                x4[:end],
                x5[:end],
                x8[:end],
                x9[:end],
            ),
        )

    results = {}

    for i, ref_id in enumerate(x2[:end]):
        if x1[i] < 0.01 and x10[i] <= 0 and x11[i] <= 0:
            pass
        else:
            results[ref_id] = {
                "final": {
                    "pi": round(x1[i], 10),
                    "best": round(x6[i], 10),
                    "high": round(x10[i], 10),
                    "low": round(x11[i], 10),
                    "reads": int(x7[i]),
                },
                "initial": {
                    "pi": round(x3[i], 10),
                    "best": round(x4[i], 10),
                    "high": round(x8[i], 10),
                    "low": round(x9[i], 10),
                    "reads": int(x5[i]),
                },
            }

    return results


def run_pathoscope(
    bam_path: Path,
    p_score_cutoff: float,
):
    """Run Pathoscope on an alignment file.

    Returns PathoscopeResults containing EM results and coverage data.

    :param alignment_path: The path to the SAM or BAM file.
    :param p_score_cutoff: The minimum allowed ``p_score`` for an alignment.
    """
    return run_expectation_maximization(
        str(bam_path),
        p_score_cutoff,
    )
