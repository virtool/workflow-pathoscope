import asyncio
import csv
import gzip
import json
import re
from pathlib import Path
from tempfile import TemporaryDirectory

from virtool.caches.utils import derive_key
from virtool.workflow.data.cache import CacheHit, WorkflowCache
from virtool.workflow.errors import SubprocessFailedError
from virtool.workflow.runtime.run_subprocess import RunSubprocess
from virtool.workflow.utils import get_workflow_version

from workflow_pathoscope.rust import PathoscopeResults, run_expectation_maximization


BOWTIE2_BUILD_TOOL = "bowtie2-build"
CD_HIT_EST_TOOL = "cd-hit-est"
CD_HIT_EST_IDENTITY = "0.99"
WORKFLOW_NAME = "pathoscope"
WORKFLOW_VERSION = get_workflow_version()


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


async def get_cd_hit_est_version(run_subprocess: RunSubprocess) -> str:
    output = []

    async def collect_output(line: bytes) -> None:
        output.append(line.decode())

    try:
        await run_subprocess(
            [CD_HIT_EST_TOOL, "-h"],
            stderr_handler=collect_output,
            stdout_handler=collect_output,
        )
    except SubprocessFailedError:
        # cd-hit-est -h prints help/version text and exits with code 1.
        pass

    match = re.search(r"\bCD-HIT\s+version\s+([^\s]+)", "".join(output))

    if match is None:
        raise ValueError("Could not parse cd-hit-est version")

    return match.group(1)


async def get_reference_collapse_cache_params(
    parent_id: str,
    run_subprocess: RunSubprocess,
) -> dict[str, str]:
    tool_version = await get_cd_hit_est_version(run_subprocess)

    return {
        "index_kind": "collapsed_reference",
        "workflow": WORKFLOW_NAME,
        "workflow_version": WORKFLOW_VERSION,
        "parent_id": parent_id,
        "source": "index_json",
        "tool_name": CD_HIT_EST_TOOL,
        "tool_version": tool_version,
        "identity": CD_HIT_EST_IDENTITY,
    }


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


def _get_otu_schema_segment_names(otu: dict) -> set[str]:
    return {str(segment["name"]) for segment in otu["schema"]}


def _get_schema_sequence_segment_key(
    otu: dict,
    sequence: dict,
    valid_schema_segments: set[str],
) -> str:
    segment_key = sequence["segment"]

    if segment_key not in valid_schema_segments:
        raise ValueError(
            f"Sequence {sequence['_id']} uses segment {segment_key!r}, which is not "
            f"defined in OTU {otu['_id']} schema"
        )

    return segment_key


def _write_fasta(sequences: list[dict], path: Path) -> None:
    with open(path, "w") as handle:
        for sequence in sequences:
            handle.write(f">{sequence['_id']}\n{sequence['sequence']}\n")


def _parse_cd_hit_clusters(cluster_path: Path) -> dict[str, str]:
    representatives_by_sequence_id = {}
    cluster_sequence_ids = []
    representative_id = None

    def flush_cluster() -> None:
        if representative_id is None:
            return

        for sequence_id in cluster_sequence_ids:
            representatives_by_sequence_id[sequence_id] = representative_id

    with open(cluster_path) as handle:
        for line in handle:
            line = line.strip()

            if line.startswith(">Cluster"):
                flush_cluster()
                cluster_sequence_ids = []
                representative_id = None
                continue

            match = re.search(r">(.+?)\.\.\.", line)

            if match is None:
                continue

            sequence_id = match.group(1)
            cluster_sequence_ids.append(sequence_id)

            if line.endswith("*"):
                representative_id = sequence_id

    flush_cluster()

    return representatives_by_sequence_id


def _build_representative_set(
    isolate: dict,
    representatives_by_sequence_id: dict[str, str],
) -> frozenset[str]:
    return frozenset(
        representatives_by_sequence_id[sequence["_id"]]
        for sequence in isolate["sequences"]
    )


async def _collapse_reference_segment(
    segment_input_path: Path,
    segment_output_path: Path,
    sequences: list[dict],
    run_subprocess: RunSubprocess,
) -> dict[str, str]:
    await asyncio.to_thread(_write_fasta, sequences, segment_input_path)

    await run_subprocess(
        [
            CD_HIT_EST_TOOL,
            "-i",
            str(segment_input_path),
            "-o",
            str(segment_output_path),
            "-c",
            CD_HIT_EST_IDENTITY,
            "-T",
            "1",
            "-M",
            "0",
            "-d",
            "0",
        ],
    )

    return _parse_cd_hit_clusters(segment_output_path.with_suffix(".cdhit.clstr"))


async def collapse_reference_json(
    json_path: Path,
    target_path: Path,
    proc: int,
    run_subprocess: RunSubprocess,
) -> dict[str, int]:
    """Collapse redundant isolates in a reference JSON using cd-hit-est clusters."""
    with _open_json(json_path) as handle:
        reference_data = json.load(handle)

    otus = _get_reference_otus(reference_data)
    before_count = 0
    after_count = 0

    with TemporaryDirectory(prefix="pathoscope-collapse-") as temp_dir:
        temp_path = Path(temp_dir)
        semaphore = asyncio.Semaphore(proc)

        async def collapse_segment(
            segment_input_path: Path,
            segment_output_path: Path,
            sequences: list[dict],
        ) -> dict[str, str]:
            async with semaphore:
                return await _collapse_reference_segment(
                    segment_input_path,
                    segment_output_path,
                    sequences,
                    run_subprocess,
                )

        for otu in otus:
            sequences_by_segment = {}
            valid_schema_segments = _get_otu_schema_segment_names(otu)
            before_count += len(otu["isolates"])
            for isolate in otu["isolates"]:
                for sequence in isolate["sequences"]:
                    segment_key = _get_schema_sequence_segment_key(
                        otu,
                        sequence,
                        valid_schema_segments,
                    )

                    sequences_by_segment.setdefault(segment_key, []).append(sequence)

            representatives_by_sequence_id = {}

            sorted_segment_sequences = sorted(
                sequences_by_segment.items(),
                key=lambda item: item[0],
            )
            segment_tasks = []

            for segment_name, sequences in sorted_segment_sequences:
                segment_input_path = (
                    temp_path / f"otu-{otu['_id']}-segment-{segment_name}.fa"
                )
                segment_output_path = (
                    temp_path / f"otu-{otu['_id']}-segment-{segment_name}.cdhit"
                )

                segment_tasks.append(
                    collapse_segment(
                        segment_input_path,
                        segment_output_path,
                        sequences,
                    )
                )

            for representatives in await asyncio.gather(*segment_tasks):
                representatives_by_sequence_id.update(representatives)

            default_sequence_ids = {
                sequence["_id"]
                for isolate in otu["isolates"]
                if isolate["default"]
                for sequence in isolate["sequences"]
            }
            seen_representative_sets = set()
            collapsed_isolates = []

            for isolate in otu["isolates"]:
                representative_set = _build_representative_set(
                    isolate,
                    representatives_by_sequence_id,
                )

                first_for_set = representative_set not in seen_representative_sets
                seen_representative_sets.add(representative_set)

                contains_default_sequence = any(
                    sequence["_id"] in default_sequence_ids
                    for sequence in isolate["sequences"]
                )

                if isolate["default"] or contains_default_sequence or first_for_set:
                    collapsed_isolates.append(isolate)

            otu["isolates"] = collapsed_isolates
            after_count += len(collapsed_isolates)

    await asyncio.to_thread(target_path.parent.mkdir, parents=True, exist_ok=True)

    with open(target_path, "w") as handle:
        json.dump(reference_data, handle)

    return {
        "isolate_count_before": before_count,
        "isolate_count_after": after_count,
        "isolate_count_removed": before_count - after_count,
    }


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
