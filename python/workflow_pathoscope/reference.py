import asyncio
import re
from collections import Counter
from collections.abc import Awaitable, Callable
from pathlib import Path
from tempfile import TemporaryDirectory

from virtool.workflow.data.indexes import WFIndex
from virtool.workflow.errors import SubprocessFailedError
from virtool.workflow.runtime.run_subprocess import RunSubprocess
from virtool.workflow.utils import get_workflow_version


CD_HIT_EST_TOOL = "cd-hit-est"
CD_HIT_EST_IDENTITY = "0.99"
WORKFLOW_NAME = "pathoscope"
WORKFLOW_VERSION = get_workflow_version()


async def get_cd_hit_est_version(run_subprocess: RunSubprocess) -> str:
    output: list[str] = []

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
        "source": "index_sqlite",
        "tool_name": CD_HIT_EST_TOOL,
        "tool_version": tool_version,
        "identity": CD_HIT_EST_IDENTITY,
    }


def _validate_isolate_sequence_segments_match_schema(
    otu: dict,
    isolate: dict,
    schema_segments: set[str],
) -> None:
    sequence_segments = [sequence["segment"] for sequence in isolate["sequences"]]
    unknown_segments = set(sequence_segments) - schema_segments

    if unknown_segments:
        raise ValueError(
            f"Isolate {isolate['id']} sequence segments "
            f"{sorted(unknown_segments, key=lambda segment: (segment is not None, segment or ''))!r} "
            f"are not defined in OTU {otu['id']} schema "
            f"segments {sorted(schema_segments)!r}"
        )

    duplicate_segments = {
        segment for segment, count in Counter(sequence_segments).items() if count > 1
    }

    if duplicate_segments:
        raise ValueError(
            f"Isolate {isolate['id']} has multiple sequences for OTU {otu['id']} "
            f"schema segments {sorted(duplicate_segments)!r}"
        )


def _validate_unsegmented_isolate(otu: dict, isolate: dict) -> None:
    sequence_count = len(isolate["sequences"])

    if sequence_count != 1:
        raise ValueError(
            f"Schema-less OTU {otu['id']} has multiple isolates, so isolate "
            f"{isolate['id']} must contain exactly one sequence; found {sequence_count}"
        )


def _write_fasta(sequences: list[dict], path: Path) -> None:
    with open(path, "w") as handle:
        for sequence in sequences:
            handle.write(f">{sequence['id']}\n{sequence['sequence']}\n")


def _parse_cd_hit_clusters(cluster_path: Path) -> dict[str, str]:
    representatives_by_sequence_id: dict[str, str] = {}
    cluster_sequence_ids: list[str] = []
    representative_id: str | None = None

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
        representatives_by_sequence_id[sequence["id"]]
        for sequence in isolate["sequences"]
    )


def _validate_and_group_otu_sequences_by_segment(otu: dict) -> dict[str, list[dict]]:
    sequences_by_segment: dict[str, list[dict]] = {}
    schema_segments = {segment["name"] for segment in otu["schema"]}

    for isolate in otu["isolates"]:
        if schema_segments:
            _validate_isolate_sequence_segments_match_schema(
                otu,
                isolate,
                schema_segments,
            )
            for sequence in isolate["sequences"]:
                sequences_by_segment.setdefault(sequence["segment"], []).append(
                    sequence,
                )
        else:
            _validate_unsegmented_isolate(otu, isolate)
            sequences_by_segment.setdefault("", []).extend(isolate["sequences"])

    return sequences_by_segment


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


async def _collapse_otu_segments(
    otu: dict,
    temp_path: Path,
    collapse_segment: Callable[[Path, Path, list[dict]], Awaitable[dict[str, str]]],
) -> dict[str, str]:
    representatives_by_sequence_id: dict[str, str] = {}
    sequences_by_segment = _validate_and_group_otu_sequences_by_segment(otu)

    segment_tasks: list[Awaitable[dict[str, str]]] = []
    sorted_segment_sequences = sorted(
        sequences_by_segment.items(),
        key=lambda item: item[0],
    )

    for segment_name, sequences in sorted_segment_sequences:
        otu_id = otu["id"]
        segment_input_path = temp_path / f"otu-{otu_id}-segment-{segment_name}.fa"
        segment_output_path = temp_path / f"otu-{otu_id}-segment-{segment_name}.cdhit"

        segment_tasks.append(
            collapse_segment(
                segment_input_path,
                segment_output_path,
                sequences,
            )
        )

    for representatives in await asyncio.gather(*segment_tasks):
        representatives_by_sequence_id.update(representatives)

    return representatives_by_sequence_id


async def _collapse_otu_reference(
    otu: dict,
    temp_path: Path,
    collapse_segment: Callable[[Path, Path, list[dict]], Awaitable[dict[str, str]]],
) -> dict:
    if not otu["schema"] and len(otu["isolates"]) == 1:
        return otu

    representatives_by_sequence_id = await _collapse_otu_segments(
        otu,
        temp_path,
        collapse_segment,
    )

    default_sequence_ids = {
        sequence["id"]
        for isolate in otu["isolates"]
        if isolate["default"]
        for sequence in isolate["sequences"]
    }
    seen_representative_sets: set[frozenset[str]] = set()
    collapsed_isolates: list[dict] = []

    for isolate in otu["isolates"]:
        representative_set = _build_representative_set(
            isolate,
            representatives_by_sequence_id,
        )

        first_for_set = representative_set not in seen_representative_sets
        seen_representative_sets.add(representative_set)

        contains_default_sequence = any(
            sequence["id"] in default_sequence_ids for sequence in isolate["sequences"]
        )

        if isolate["default"] or contains_default_sequence or first_for_set:
            collapsed_isolates.append(isolate)

    return {**otu, "isolates": collapsed_isolates}


async def collapse_reference_index(
    index: WFIndex,
    target_path: Path,
    proc: int,
    run_subprocess: RunSubprocess,
) -> dict[str, int]:
    """Collapse redundant isolates into a new workflow SQLite index."""
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

        collapsed_otus: list[dict] = []

        async for otu in index.iter_otus():
            collapsed_otu = await _collapse_otu_reference(
                otu,
                temp_path,
                collapse_segment,
            )
            collapsed_otus.append(collapsed_otu)
            before_count += len(otu["isolates"])
            after_count += len(collapsed_otu["isolates"])

    try:
        reference = await index.get_reference_metadata()
    except ValueError:
        reference = None

    await WFIndex.create(index.id, target_path, reference, collapsed_otus)

    return {
        "isolate_count_before": before_count,
        "isolate_count_after": after_count,
        "isolate_count_removed": before_count - after_count,
    }
