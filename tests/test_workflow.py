import asyncio
import re
import shutil
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest
from structlog import get_logger
from syrupy import SnapshotAssertion
from virtool.caches.utils import derive_key
from virtool.workflow import RunSubprocess
from virtool.workflow.data.analyses import WFAnalysis
from virtool.workflow.data.cache import WorkflowCache
from virtool.workflow.data.indexes import WFIndex
from virtool.workflow.data.samples import WFSample
from virtool.workflow.data.subtractions import WFSubtraction
from virtool.workflow.errors import JobsAPINotFoundError
from virtool.workflow.pytest_plugin import WorkflowData
from workflow_pathoscope.reference import (
    CD_HIT_EST_IDENTITY,
    _prepare_otu_collapse,
    collapse_reference_index,
    get_reference_collapse_cache_params,
)
from workflow_pathoscope.utils import get_mapping_index_cache_params

from fixtures import (
    get_collapsed_reference_path,
    get_isolate_bam_path,
    get_isolate_fastq_path,
    get_isolate_index_path,
    get_isolate_path,
    get_reference_index_path,
    get_subtracted_bam_path,
    get_subtraction_indexes_path,
)
from workflow import (
    build_isolate_index,
    collapse_reference,
    create_reference_index,
    create_subtraction_index,
    eliminate_subtraction,
    get_subtraction_index_path,
    map_default_isolates,
    map_isolates,
    reassignment,
)

BOWTIE2_INDEX_SUFFIXES = (
    "1.bt2",
    "2.bt2",
    "3.bt2",
    "4.bt2",
    "rev.1.bt2",
    "rev.2.bt2",
)
TOOL_VERSION_PATTERN = re.compile(r"\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?")


@pytest.fixture()
def work_path(tmpdir):
    path = Path(tmpdir) / "work"
    path.mkdir(parents=True)

    return path


@pytest.fixture()
def reference_index_path(work_path: Path) -> Path:
    return get_reference_index_path(work_path)


@pytest.fixture()
def collapsed_reference_path(work_path: Path) -> Path:
    return get_collapsed_reference_path(work_path)


@pytest.fixture()
def subtraction_indexes_path(work_path: Path) -> Path:
    return get_subtraction_indexes_path(work_path)


@pytest.fixture()
def subtraction_index_path(
    subtractions: list[WFSubtraction],
    subtraction_indexes_path: Path,
) -> Path:
    return get_subtraction_index_path(subtraction_indexes_path, subtractions[0].id)


@pytest.fixture()
def isolate_path(work_path: Path) -> Path:
    path = get_isolate_path(work_path)
    path.mkdir()

    return path


@pytest.fixture()
def isolate_fastq_path(isolate_path: Path) -> Path:
    return get_isolate_fastq_path(isolate_path)


@pytest.fixture()
def isolate_index_path(isolate_path: Path) -> Path:
    return get_isolate_index_path(isolate_path)


@pytest.fixture()
def isolate_bam_path(isolate_path: Path) -> Path:
    return get_isolate_bam_path(isolate_path)


@pytest.fixture()
def subtracted_bam_path(work_path: Path) -> Path:
    return get_subtracted_bam_path(work_path)


class _FakeWorkflowCacheAPI:
    """Fake only the API calls used by the real workflow cache."""

    def __init__(
        self,
        work_dir: Path,
        *,
        put_exception: Exception | None = None,
        put_created: bool | None = None,
    ) -> None:
        self.work_dir = work_dir
        self.put_exception = put_exception
        self.put_created = put_created
        self.stored: dict[str, tuple[Path, dict | None]] = {}

    async def get_cache(self, key: str, dest: Path) -> None:
        try:
            source, _ = self.stored[key]
        except KeyError:
            raise JobsAPINotFoundError from None

        dest.parent.mkdir(parents=True, exist_ok=True)
        await asyncio.to_thread(shutil.copyfile, source, dest)

    async def put_cache(
        self,
        key: str,
        path: Path,
        params: dict | None = None,
    ) -> bool:
        if self.put_exception:
            raise self.put_exception

        if key in self.stored:
            return False

        stored_path = self.work_dir / "stored" / key / "cache.tar"
        stored_path.parent.mkdir(parents=True, exist_ok=True)
        await asyncio.to_thread(shutil.copyfile, path, stored_path)
        self.stored[key] = (stored_path, params)

        if self.put_created is not None:
            return self.put_created

        return True


@pytest.fixture()
def workflow_cache(tmp_path: Path) -> WorkflowCache:
    api = _FakeWorkflowCacheAPI(tmp_path / "fake_workflow_cache_api")
    return WorkflowCache(api, tmp_path / "workflow_cache")


def read_directory_bytes(path: Path) -> dict[str, bytes]:
    return {child.name: child.read_bytes() for child in sorted(path.iterdir())}


def bowtie2_bundle_bytes(prefix: str, content: bytes | None = None) -> dict[str, bytes]:
    return {
        f"{prefix}.{suffix}": content or f"{prefix}.{suffix}".encode()
        for suffix in BOWTIE2_INDEX_SUFFIXES
    }


def write_bowtie2_bundle(
    path: Path,
    prefix: str,
    content: bytes = b"cached",
    extra_files: dict[str, bytes] | None = None,
):
    path.mkdir(parents=True)

    for name, file_content in bowtie2_bundle_bytes(prefix, content).items():
        (path / name).write_bytes(file_content)

    for name, file_content in (extra_files or {}).items():
        (path / name).write_bytes(file_content)


def assert_bowtie2_index_exists(prefix: Path):
    assert read_directory_bytes(prefix.parent).keys() == {
        f"{prefix.name}.{suffix}" for suffix in BOWTIE2_INDEX_SUFFIXES
    }

    for suffix in BOWTIE2_INDEX_SUFFIXES:
        assert (prefix.parent / f"{prefix.name}.{suffix}").stat().st_size > 0


def assert_cache_params(
    params: dict[str, str],
    expected: dict[str, str],
) -> None:
    assert params.keys() == {*expected.keys(), "tool_version"}
    assert {
        key: value for key, value in params.items() if key != "tool_version"
    } == expected
    assert TOOL_VERSION_PATTERN.fullmatch(params["tool_version"])


async def as_async_iterable(items):
    for item in items:
        yield item


async def create_wf_index(id_, path, reference, otus):
    return await WFIndex.create(id_, path, reference, as_async_iterable(otus))


async def create_test_index(path: Path) -> WFIndex:
    def sequence(id_: str, value: str) -> dict:
        return {
            "id": id_,
            "accession": id_,
            "definition": id_,
            "host": "",
            "segment": None,
            "sequence": value,
        }

    def isolate(id_: str, default: bool, sequences: list[dict]) -> dict:
        return {
            "id": id_,
            "default": default,
            "sequences": sequences,
            "source_name": id_,
            "source_type": "isolate",
        }

    return await create_wf_index(
        "test-index",
        path,
        None,
        [
            {
                "id": "default-otu",
                "abbreviation": "DEFAULT",
                "isolates": [
                    isolate(
                        "default",
                        True,
                        [
                            sequence("default-a", "ACGT"),
                            sequence("default-b", "TTA"),
                        ],
                    ),
                    isolate(
                        "non-default",
                        False,
                        [sequence("non-default", "GGGG")],
                    ),
                ],
                "name": "Default OTU",
                "schema": [],
                "taxid": None,
                "version": 1,
            },
            {
                "id": "non-default-otu",
                "abbreviation": "NONDEFAULT",
                "isolates": [
                    isolate(
                        "non-default-only",
                        False,
                        [sequence("non-default-only", "CCCC")],
                    ),
                ],
                "name": "Non-default OTU",
                "schema": [],
                "taxid": None,
                "version": 1,
            },
        ],
    )


@pytest.fixture()
def analysis(workflow_data: WorkflowData, mocker):
    analysis_ = mocker.Mock(WFAnalysis)

    analysis_.id = workflow_data.analysis.id
    analysis_.workflow = "pathoscope_bowtie"
    analysis_.ready = False
    analysis_.sample = workflow_data.analysis.sample

    return analysis_


@pytest.fixture()
async def index(workflow_data: WorkflowData, work_path: Path):
    sequence_otu_map = {
        "NC_016509": "foobar",
        "NC_001948": "foobar",
        "13TF149_Reovirus_TF1_Seg06": "reo",
        "13TF149_Reovirus_TF1_Seg03": "reo",
        "13TF149_Reovirus_TF1_Seg07": "reo",
        "13TF149_Reovirus_TF1_Seg02": "reo",
        "13TF149_Reovirus_TF1_Seg08": "reo",
        "13TF149_Reovirus_TF1_Seg11": "reo",
        "13TF149_Reovirus_TF1_Seg04": "reo",
        "NC_004667": "foobar",
        "NC_003347": "foobar",
        "NC_003615": "foobar",
        "NC_003689": "foobar",
        "NC_011552": "foobar",
        "KX109927": "baz",
        "NC_008039": "foobar",
        "NC_015782": "foobar",
        "NC_016416": "foobar",
        "NC_003623": "foobar",
        "NC_008038": "foobar",
        "NC_001836": "foobar",
        "JQ080272": "baz",
        "NC_017938": "foobar",
        "NC_008037": "foobar",
        "NC_007448": "foobar",
    }
    manifest = {"foobar": 10, "reo": 5, "baz": 6}
    otus = []

    for otu_id, version in manifest.items():
        otus.append(
            {
                "id": otu_id,
                "abbreviation": otu_id.upper(),
                "isolates": [
                    {
                        "default": True,
                        "id": f"{otu_id}-isolate",
                        "sequences": [
                            {
                                "id": sequence_id,
                                "accession": sequence_id,
                                "definition": sequence_id,
                                "host": "",
                                "segment": None,
                                "sequence": "ACGT",
                            }
                            for sequence_id, sequence_otu_id in sequence_otu_map.items()
                            if sequence_otu_id == otu_id
                        ],
                        "source_name": otu_id,
                        "source_type": "isolate",
                    },
                ],
                "name": otu_id,
                "schema": [],
                "taxid": None,
                "version": version,
            },
        )

    return await create_wf_index(
        workflow_data.index.id,
        work_path / "indexes" / str(workflow_data.index.id) / "index.sqlite",
        None,
        otus,
    )


def make_index_sequence(
    sequence_id: str,
    *,
    segment: object,
    sequence: str = "ACGT",
) -> dict:
    return {
        "id": sequence_id,
        "accession": sequence_id,
        "definition": sequence_id,
        "host": "",
        "segment": segment,
        "sequence": sequence,
    }


def make_index_isolate(
    isolate_id: str,
    *,
    sequences: list[dict],
    default: bool = False,
) -> dict:
    return {
        "id": isolate_id,
        "default": default,
        "sequences": sequences,
        "source_name": isolate_id,
        "source_type": "isolate",
    }


def make_index_otu(
    otu_id: str,
    *,
    isolates: list[dict],
    schema: list[str],
) -> dict:
    return {
        "id": otu_id,
        "abbreviation": "COLLAPSE",
        "name": "Collapse OTU",
        "taxid": None,
        "version": 1,
        "schema": [{"name": segment} for segment in schema],
        "isolates": isolates,
    }


def make_redundant_index_otu() -> dict:
    sequence_data = {
        "default": (
            "ACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCA",
        ),
        "representative-1": (
            "ACGTACGTACGTACGTACGT",
            "TGCAGGATCGTTTAACGTAG",
        ),
        "representative-2": (
            "CCCCAAAAGGGGTTTTCCCC",
            "AAAACCCCGGGGTTTTAAAA",
        ),
        "unique-combo": (
            "CCCCAAAAGGGGTTTTCCCC",
            "TGCAGGATCGTTTAACGTAG",
        ),
        "duplicate": (
            "CCCCAAAAGGGGTTTTCCCC",
            "AAAACCCCGGGGTTTTAAAA",
        ),
    }

    return make_index_otu(
        "collapse-otu",
        schema=["a", "b"],
        isolates=[
            make_index_isolate(
                isolate_id,
                default=isolate_id == "default",
                sequences=[
                    make_index_sequence(
                        f"{isolate_id}-a",
                        segment="a",
                        sequence=sequences[0],
                    ),
                    make_index_sequence(
                        f"{isolate_id}-b",
                        segment="b",
                        sequence=sequences[1],
                    ),
                ],
            )
            for isolate_id, sequences in sequence_data.items()
        ],
    )


@pytest.fixture()
async def redundant_index(workflow_data: WorkflowData, tmp_path: Path) -> WFIndex:
    return await create_wf_index(
        workflow_data.index.id,
        tmp_path / "redundant-index.sqlite",
        None,
        [make_redundant_index_otu()],
    )


async def test_collapse_reference_hit(
    collapsed_reference_path: Path,
    index: WFIndex,
    run_subprocess: RunSubprocess,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
):
    source = tmp_path / collapsed_reference_path.parent.name
    source.mkdir()
    await create_test_index(source / collapsed_reference_path.name)
    params = await get_reference_collapse_cache_params(index.id, run_subprocess)
    key = derive_key(params)
    assert await workflow_cache.put(key, source, params)

    logger = get_logger("test")

    result_path = await collapse_reference(
        workflow_cache,
        collapsed_reference_path,
        index,
        logger,
        4,
        run_subprocess,
    )

    assert result_path == collapsed_reference_path
    assert read_directory_bytes(
        collapsed_reference_path.parent
    ) == read_directory_bytes(source)


async def test_collapse_reference_miss_retains_required_isolates(
    collapsed_reference_path: Path,
    redundant_index: WFIndex,
    run_subprocess: RunSubprocess,
    workflow_cache: WorkflowCache,
):
    logger = get_logger("test")

    result_path = await collapse_reference(
        workflow_cache,
        collapsed_reference_path,
        redundant_index,
        logger,
        4,
        run_subprocess,
    )

    params = await get_reference_collapse_cache_params(
        redundant_index.id,
        run_subprocess,
    )

    assert result_path == collapsed_reference_path
    assert_cache_params(
        params,
        {
            "index_kind": "collapsed_reference",
            "workflow": "pathoscope",
            "workflow_version": "UNKNOWN",
            "parent_id": redundant_index.id,
            "source": "index_sqlite",
            "tool_name": "cd-hit-est",
            "identity": "0.99",
        },
    )

    collapsed_index = WFIndex.load(redundant_index.id, collapsed_reference_path)
    collapsed_otus = [otu async for otu in collapsed_index.iter_otus()]

    assert [isolate["id"] for isolate in collapsed_otus[0]["isolates"]] == [
        "default",
        "representative-1",
        "representative-2",
        "unique-combo",
    ]
    assert [
        (isolate["sequences"][0]["id"], isolate["sequences"][1]["id"])
        for isolate in collapsed_otus[0]["isolates"]
    ] == [
        ("default-a", "default-b"),
        ("representative-1-a", "representative-1-b"),
        ("representative-2-a", "representative-2-b"),
        ("unique-combo-a", "unique-combo-b"),
    ]
    assert not (collapsed_reference_path.parent / "collapse-manifest.json").exists()


async def test_collapse_reference_index_outputs_collapsed_reference_fasta(
    redundant_index: WFIndex,
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"
    default_fasta_path = tmp_path / "default.fa"
    isolate_fasta_path = tmp_path / "isolates.fa"

    assert await collapse_reference_index(
        redundant_index,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 4,
        "isolate_count_removed": 1,
        "otu_count_collapsed": 1,
        "otu_count_unchanged": 0,
        "otu_count_skipped": 0,
    }

    collapsed_index = WFIndex.load(redundant_index.id, collapsed_path)

    await collapsed_index.write_fasta(
        default_fasta_path,
        collapsed_index.iter_default_sequences(),
    )
    await collapsed_index.write_fasta(
        isolate_fasta_path,
        collapsed_index.iter_otu_sequences("collapse-otu"),
    )

    assert default_fasta_path.read_text() == (
        ">default-a\nACGTACGTACGTACGTACGT\n>default-b\nTGCATGCATGCATGCATGCA\n"
    )
    assert "duplicate-a" not in isolate_fasta_path.read_text()
    assert "duplicate-b" not in isolate_fasta_path.read_text()


async def test_collapse_reference_index_allows_isolates_missing_schema_segments(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    otu = make_redundant_index_otu()
    default_isolate = next(
        isolate for isolate in otu["isolates"] if isolate["id"] == "default"
    )
    default_isolate["sequences"] = [
        sequence
        for sequence in default_isolate["sequences"]
        if sequence["segment"] == "b"
    ]
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    assert await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 4,
        "isolate_count_removed": 1,
        "otu_count_collapsed": 1,
        "otu_count_unchanged": 0,
        "otu_count_skipped": 0,
    }

    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    collapsed_otus = [otu async for otu in collapsed_index.iter_otus()]

    assert [
        sequence["id"] for sequence in collapsed_otus[0]["isolates"][0]["sequences"]
    ] == ["default-b"]


async def test_collapse_reference_index_allows_unsegmented_isolates(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    otu = make_index_otu(
        "collapse-otu",
        schema=[],
        isolates=[
            make_index_isolate(
                "default",
                default=True,
                sequences=[
                    make_index_sequence(
                        "default",
                        segment=None,
                        sequence="ACGTACGTACGTACGTACGT",
                    ),
                ],
            ),
            make_index_isolate(
                "representative-1",
                sequences=[
                    make_index_sequence(
                        "representative-1",
                        segment="ignored",
                        sequence="ACGTACGTACGTACGTACGT",
                    ),
                ],
            ),
            make_index_isolate(
                "representative-2",
                sequences=[
                    make_index_sequence(
                        "representative-2",
                        segment="ignored",
                        sequence="CCCCAAAAGGGGTTTTCCCC",
                    ),
                ],
            ),
            make_index_isolate(
                "unique-combo",
                sequences=[
                    make_index_sequence(
                        "unique-combo",
                        segment="ignored",
                        sequence="CCCCAAAAGGGGTTTTCCCC",
                    ),
                ],
            ),
            make_index_isolate(
                "duplicate",
                sequences=[
                    make_index_sequence(
                        "duplicate",
                        segment="ignored",
                        sequence="CCCCAAAAGGGGTTTTCCCC",
                    ),
                ],
            ),
        ],
    )

    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    assert await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 2,
        "isolate_count_removed": 3,
        "otu_count_collapsed": 1,
        "otu_count_unchanged": 0,
        "otu_count_skipped": 0,
    }

    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    collapsed_otus = [otu async for otu in collapsed_index.iter_otus()]

    assert collapsed_otus[0]["schema"] == []
    assert [isolate["id"] for isolate in collapsed_otus[0]["isolates"]] == [
        "default",
        "representative-2",
    ]


async def test_collapse_reference_index_rejects_multiple_sequences_for_unsegmented_otu(
    mocker,
    tmp_path: Path,
):
    otu = make_index_otu(
        "unsegmented",
        schema=[],
        isolates=[
            make_index_isolate(
                "multiple-sequences",
                default=True,
                sequences=[
                    make_index_sequence("sequence-1", segment=None),
                    make_index_sequence("sequence-2", segment=None),
                ],
            ),
            make_index_isolate(
                "single-sequence",
                sequences=[make_index_sequence("sequence-3", segment=None)],
            ),
        ],
    )

    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    run_subprocess = mocker.AsyncMock()
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    assert await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 2,
        "isolate_count_after": 2,
        "isolate_count_removed": 0,
        "otu_count_collapsed": 0,
        "otu_count_unchanged": 0,
        "otu_count_skipped": 1,
    }

    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    assert [otu async for otu in collapsed_index.iter_otus()] == [
        otu async for otu in source_index.iter_otus()
    ]
    run_subprocess.assert_not_awaited()


async def test_collapse_reference_index_preserves_single_unsegmented_isolate(
    mocker,
    tmp_path: Path,
):
    otu = make_index_otu(
        "unsegmented",
        schema=[],
        isolates=[
            make_index_isolate(
                "only-isolate",
                default=True,
                sequences=[
                    make_index_sequence("sequence-1", segment=None),
                    make_index_sequence("sequence-2", segment=None),
                ],
            ),
        ],
    )
    run_subprocess = mocker.AsyncMock()
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    assert await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 1,
        "isolate_count_after": 1,
        "isolate_count_removed": 0,
        "otu_count_collapsed": 0,
        "otu_count_unchanged": 1,
        "otu_count_skipped": 0,
    }

    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    collapsed_otus = [otu async for otu in collapsed_index.iter_otus()]

    assert [
        sequence["id"] for sequence in collapsed_otus[0]["isolates"][0]["sequences"]
    ] == ["sequence-1", "sequence-2"]
    run_subprocess.assert_not_awaited()


async def test_collapse_reference_index_preserves_single_schema_isolate_without_validation(
    mocker,
    tmp_path: Path,
):
    otu = make_index_otu(
        "segmented",
        schema=["a", "b"],
        isolates=[
            make_index_isolate(
                "only-isolate",
                default=True,
                sequences=[
                    make_index_sequence("sequence-a", segment=None),
                    make_index_sequence("sequence-b", segment="b"),
                ],
            ),
        ],
    )
    run_subprocess = mocker.AsyncMock()
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    assert await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 1,
        "isolate_count_after": 1,
        "isolate_count_removed": 0,
        "otu_count_collapsed": 0,
        "otu_count_unchanged": 1,
        "otu_count_skipped": 0,
    }

    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    assert [otu async for otu in collapsed_index.iter_otus()] == [
        otu async for otu in source_index.iter_otus()
    ]
    run_subprocess.assert_not_awaited()


@pytest.mark.parametrize("segment", [None, ""])
async def test_collapse_reference_index_skips_invalid_schema_segment_names(
    mocker,
    segment,
    tmp_path: Path,
):
    otu = make_index_otu(
        "segmented",
        schema=["a", "b"],
        isolates=[
            make_index_isolate(
                "invalid",
                default=True,
                sequences=[
                    make_index_sequence("invalid-a", segment=segment),
                    make_index_sequence("invalid-b", segment="b"),
                ],
            ),
            make_index_isolate(
                "valid",
                sequences=[
                    make_index_sequence("valid-a", segment="a"),
                    make_index_sequence("valid-b", segment="b"),
                ],
            ),
        ],
    )
    run_subprocess = mocker.AsyncMock()
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    summary = await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    )

    assert summary["otu_count_skipped"] == 1
    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    assert [otu async for otu in collapsed_index.iter_otus()] == [
        otu async for otu in source_index.iter_otus()
    ]
    run_subprocess.assert_not_awaited()


def test_otu_collapse_preparation_rejects_non_string_schema_segment():
    otu = make_index_otu(
        "segmented",
        schema=["a"],
        isolates=[
            make_index_isolate(
                "invalid",
                sequences=[make_index_sequence("invalid-a", segment=42)],
            ),
        ],
    )

    preparation = _prepare_otu_collapse(otu)

    assert preparation.eligible is False
    assert preparation.sequences_by_segment == {}


def test_otu_collapse_preparation_rejects_empty_schemaless_isolate():
    otu = make_index_otu(
        "unsegmented",
        schema=[],
        isolates=[
            make_index_isolate("empty", sequences=[]),
            make_index_isolate(
                "valid",
                sequences=[make_index_sequence("valid", segment=None)],
            ),
        ],
    )

    preparation = _prepare_otu_collapse(otu)

    assert preparation.eligible is False
    assert preparation.sequences_by_segment == {}


async def test_collapse_reference_index_rejects_duplicate_isolate_segments(
    mocker,
    tmp_path: Path,
):
    otu = make_index_otu(
        "segmented",
        schema=["a", "b"],
        isolates=[
            make_index_isolate(
                "duplicate",
                default=True,
                sequences=[
                    make_index_sequence("duplicate-a-1", segment="a"),
                    make_index_sequence("duplicate-a-2", segment="a"),
                ],
            ),
            make_index_isolate(
                "valid",
                sequences=[
                    make_index_sequence("valid-a", segment="a"),
                    make_index_sequence("valid-b", segment="b"),
                ],
            ),
        ],
    )
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    run_subprocess = mocker.AsyncMock()
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    summary = await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    )

    assert summary["otu_count_skipped"] == 1
    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    assert [otu async for otu in collapsed_index.iter_otus()] == [
        otu async for otu in source_index.iter_otus()
    ]
    run_subprocess.assert_not_awaited()


async def test_collapse_reference_index_rejects_unknown_isolate_segments(
    mocker,
    tmp_path: Path,
):
    otu = make_index_otu(
        "segmented",
        schema=["a", "b"],
        isolates=[
            make_index_isolate(
                "unknown",
                default=True,
                sequences=[
                    make_index_sequence("unknown-c", segment="c"),
                    make_index_sequence("unknown-b", segment="b"),
                ],
            ),
            make_index_isolate(
                "valid",
                sequences=[
                    make_index_sequence("valid-a", segment="a"),
                    make_index_sequence("valid-b", segment="b"),
                ],
            ),
        ],
    )
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [otu],
    )
    run_subprocess = mocker.AsyncMock()
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"

    summary = await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        run_subprocess,
    )

    assert summary["otu_count_skipped"] == 1
    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    assert [otu async for otu in collapsed_index.iter_otus()] == [
        otu async for otu in source_index.iter_otus()
    ]
    run_subprocess.assert_not_awaited()


async def test_collapse_reference_index_handles_mixed_otu_outcomes(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    eligible_otu = make_redundant_index_otu()
    skipped_otu = make_index_otu(
        "skipped",
        schema=["a"],
        isolates=[
            make_index_isolate(
                "skipped-invalid",
                default=True,
                sequences=[
                    make_index_sequence("skipped-invalid", segment=None),
                ],
            ),
            make_index_isolate(
                "skipped-valid",
                sequences=[
                    make_index_sequence("skipped-valid", segment="a"),
                ],
            ),
        ],
    )
    unchanged_otu = make_index_otu(
        "unchanged",
        schema=["a"],
        isolates=[
            make_index_isolate(
                "unchanged-invalid",
                default=True,
                sequences=[
                    make_index_sequence("unchanged-invalid", segment=None),
                ],
            ),
        ],
    )
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [eligible_otu, skipped_otu, unchanged_otu],
    )
    collapsed_path = tmp_path / "collapsed" / "index.sqlite"
    subprocess_commands: list[list[str]] = []

    async def tracked_run_subprocess(command, **kwargs):
        subprocess_commands.append(command)
        return await run_subprocess(command, **kwargs)

    assert await collapse_reference_index(
        source_index,
        collapsed_path,
        2,
        tracked_run_subprocess,
    ) == {
        "isolate_count_before": 8,
        "isolate_count_after": 7,
        "isolate_count_removed": 1,
        "otu_count_collapsed": 1,
        "otu_count_unchanged": 1,
        "otu_count_skipped": 1,
    }

    collapsed_index = WFIndex.load(source_index.id, collapsed_path)
    collapsed_otus = {otu["id"]: otu async for otu in collapsed_index.iter_otus()}
    source_otus = {otu["id"]: otu async for otu in source_index.iter_otus()}

    assert len(collapsed_otus["collapse-otu"]["isolates"]) == 4
    assert collapsed_otus["skipped"] == source_otus["skipped"]
    assert collapsed_otus["unchanged"] == source_otus["unchanged"]
    assert len(subprocess_commands) == 2


async def test_collapse_reference_index_propagates_subprocess_failures(
    mocker,
    tmp_path: Path,
):
    source_index = await create_wf_index(
        "test-index",
        tmp_path / "source.sqlite",
        None,
        [make_redundant_index_otu()],
    )
    run_subprocess = mocker.AsyncMock(side_effect=RuntimeError("cd-hit-est failed"))

    with pytest.raises(RuntimeError, match="cd-hit-est failed"):
        await collapse_reference_index(
            source_index,
            tmp_path / "collapsed" / "index.sqlite",
            2,
            run_subprocess,
        )


@pytest.fixture()
def sample(workflow_data: WorkflowData, example_path: Path, work_path: Path):
    workflow_data.sample.library_type = "normal"

    path = work_path / "samples" / str(workflow_data.sample.id)
    path.mkdir(parents=True)

    shutil.copyfile(example_path / "sample" / "reads_1.fq.gz", path / "reads_1.fq.gz")

    return WFSample(
        id=workflow_data.sample.id,
        library_type=workflow_data.sample.library_type,
        name=workflow_data.sample.name,
        paired=False,
        quality=workflow_data.sample.quality,
        read_paths=(path / "reads_1.fq.gz",),
    )


@pytest.fixture()
def subtractions(workflow_data: WorkflowData, example_path: Path, work_path: Path):
    subtraction_path = work_path / "subtractions" / "subtraction"
    subtraction_path.parent.mkdir(parents=True)

    shutil.copytree(
        example_path / "subtractions" / "arabidopsis_thaliana", subtraction_path
    )

    return [
        WFSubtraction(
            id=workflow_data.subtraction.id,
            files=[],
            gc=workflow_data.subtraction.gc,
            name=workflow_data.subtraction.name,
            nickname=workflow_data.subtraction.nickname,
            path=subtraction_path,
        ),
    ]


async def test_create_reference_index_hit(
    collapsed_reference_path: Path,
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
):
    await create_test_index(collapsed_reference_path)
    source = tmp_path / reference_index_path.parent.name
    write_bowtie2_bundle(
        source,
        "reference",
        b"cached-reference",
        {"cache-manifest.json": b"cached-manifest"},
    )
    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        run_subprocess,
        {
            "collapse_identity": CD_HIT_EST_IDENTITY,
            "source": "collapsed_reference",
            "selection": "default_isolates",
        },
    )
    key = derive_key(params)
    assert await workflow_cache.put(key, source, params)

    logger = get_logger("test")

    result_index_path = await create_reference_index(
        workflow_cache,
        collapsed_reference_path,
        index,
        logger,
        4,
        run_subprocess,
        reference_index_path,
    )

    assert result_index_path == reference_index_path
    assert read_directory_bytes(reference_index_path.parent) == read_directory_bytes(
        source
    )


async def test_create_reference_index_miss(
    collapsed_reference_path: Path,
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    workflow_cache: WorkflowCache,
):
    await create_test_index(collapsed_reference_path)
    logger = get_logger("test")

    result_index_path = await create_reference_index(
        workflow_cache,
        collapsed_reference_path,
        index,
        logger,
        4,
        run_subprocess,
        reference_index_path,
    )

    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        run_subprocess,
        {
            "collapse_identity": CD_HIT_EST_IDENTITY,
            "source": "collapsed_reference",
            "selection": "default_isolates",
        },
    )
    assert result_index_path == reference_index_path
    assert_cache_params(
        params,
        {
            "index_kind": "reference_mapping_index",
            "workflow": "pathoscope",
            "workflow_version": "UNKNOWN",
            "parent_id": index.id,
            "collapse_identity": "0.99",
            "source": "collapsed_reference",
            "selection": "default_isolates",
            "tool_name": "bowtie2-build",
        },
    )
    assert_bowtie2_index_exists(reference_index_path)


async def test_create_reference_index_continues_when_cache_put_is_skipped(
    collapsed_reference_path: Path,
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    await create_test_index(collapsed_reference_path)
    cache = WorkflowCache(
        _FakeWorkflowCacheAPI(
            tmp_path / "fake_workflow_cache_api",
            put_created=False,
        ),
        tmp_path / "workflow_cache",
    )
    logger = get_logger("test")

    result_index_path = await create_reference_index(
        cache,
        collapsed_reference_path,
        index,
        logger,
        4,
        run_subprocess,
        reference_index_path,
    )

    assert result_index_path == reference_index_path
    assert_bowtie2_index_exists(reference_index_path)


async def test_create_reference_index_raises_unexpected_cache_put_failure(
    collapsed_reference_path: Path,
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    await create_test_index(collapsed_reference_path)
    cache = WorkflowCache(
        _FakeWorkflowCacheAPI(
            tmp_path / "fake_workflow_cache_api",
            put_exception=RuntimeError("cache upload failed"),
        ),
        tmp_path / "workflow_cache",
    )
    logger = get_logger("test")

    with pytest.raises(RuntimeError, match="cache upload failed"):
        await create_reference_index(
            cache,
            collapsed_reference_path,
            index,
            logger,
            4,
            run_subprocess,
            reference_index_path,
        )


async def test_create_subtraction_index_hit(
    run_subprocess: RunSubprocess,
    subtractions: list[WFSubtraction],
    subtraction_index_path: Path,
    subtraction_indexes_path: Path,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
):
    source = tmp_path / subtraction_index_path.parent.name
    write_bowtie2_bundle(
        source,
        "subtraction",
        b"cached-subtraction",
        {"cache-manifest.json": b"cached-manifest"},
    )
    subtraction = subtractions[0]
    params = await get_mapping_index_cache_params(
        "subtraction_mapping_index",
        subtraction.id,
        run_subprocess,
    )
    key = derive_key(params)
    assert await workflow_cache.put(key, source, params)

    logger = get_logger("test")

    result_indexes_path = await create_subtraction_index(
        workflow_cache,
        logger,
        4,
        run_subprocess,
        subtractions,
        subtraction_indexes_path,
    )

    assert result_indexes_path == subtraction_indexes_path
    assert read_directory_bytes(subtraction_index_path.parent) == read_directory_bytes(
        source
    )


async def test_create_subtraction_index_miss(
    run_subprocess: RunSubprocess,
    subtractions: list[WFSubtraction],
    subtraction_index_path: Path,
    subtraction_indexes_path: Path,
    workflow_cache: WorkflowCache,
):
    logger = get_logger("test")
    subtraction = subtractions[0]

    result_indexes_path = await create_subtraction_index(
        workflow_cache,
        logger,
        4,
        run_subprocess,
        subtractions,
        subtraction_indexes_path,
    )

    params = await get_mapping_index_cache_params(
        "subtraction_mapping_index",
        subtraction.id,
        run_subprocess,
    )
    assert result_indexes_path == subtraction_indexes_path
    assert_cache_params(
        params,
        {
            "index_kind": "subtraction_mapping_index",
            "workflow": "pathoscope",
            "workflow_version": "UNKNOWN",
            "parent_id": subtraction.id,
            "tool_name": "bowtie2-build",
        },
    )
    assert_bowtie2_index_exists(subtraction_index_path)


async def test_map_default_isolates(
    example_path: Path,
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess,
    sample: WFSample,
    snapshot: SnapshotAssertion,
):
    reference_index_path.parent.mkdir(parents=True)

    for path in (example_path / "index").iterdir():
        if "reference" in path.name:
            shutil.copyfile(path, reference_index_path.parent / path.name)

    intermediate = SimpleNamespace(candidate_sequence_ids=set())

    logger = get_logger("test")

    await map_default_isolates(
        intermediate,
        logger,
        index,
        2,
        0.01,
        reference_index_path,
        sample,
    )

    assert sorted(intermediate.candidate_sequence_ids) == snapshot


async def test_map_isolates(
    example_path: Path,
    index: WFIndex,
    isolate_bam_path: Path,
    isolate_fastq_path: Path,
    isolate_index_path: Path,
    sample: WFSample,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
):
    for path in (example_path / "index").iterdir():
        if "reference" in path.name:
            shutil.copyfile(
                path,
                isolate_index_path.parent / path.name.replace("reference", "isolates"),
            )

    proc = 1

    intermediate = SimpleNamespace(candidate_sequence_ids={"candidate"})

    await map_isolates(
        intermediate,
        isolate_fastq_path,
        isolate_index_path,
        isolate_bam_path,
        proc,
        run_subprocess,
        sample,
    )

    # Convert BAM to text for snapshot comparison
    with pysam.AlignmentFile(isolate_bam_path, "rb") as bam:
        lines = []
        for read in bam:
            # Convert each record to SAM format text (excluding header)
            lines.append(read.to_string())

        assert sorted(lines) == snapshot

    # The test just verifies that the SAM file is written correctly


@pytest.mark.parametrize(
    "no_subtractions",
    [True, False],
    ids=["no_subtractions", "single_subtraction"],
)
async def test_eliminate_subtraction(
    example_path: Path,
    isolate_bam_path: Path,
    isolate_fastq_path: Path,
    no_subtractions: bool,
    subtracted_bam_path: Path,
    subtractions: list[WFSubtraction],
    subtraction_indexes_path: Path,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
    work_path: Path,
):
    shutil.copyfile(example_path / "to_isolates.bam", isolate_bam_path)
    shutil.copyfile(example_path / "to_isolates.fq", isolate_fastq_path)

    logger = get_logger("test")
    proc = 2
    results = {}

    if no_subtractions:
        subtractions = []

    intermediate = SimpleNamespace(candidate_sequence_ids={"candidate"})

    if subtractions:
        cached_subtraction_path = tmp_path / str(subtractions[0].id)
        shutil.copytree(subtractions[0].path, cached_subtraction_path)
        params = await get_mapping_index_cache_params(
            "subtraction_mapping_index",
            subtractions[0].id,
            run_subprocess,
        )
        key = derive_key(params)
        assert await workflow_cache.put(key, cached_subtraction_path, params)

        await create_subtraction_index(
            workflow_cache,
            logger,
            proc,
            run_subprocess,
            subtractions,
            subtraction_indexes_path,
        )

    await eliminate_subtraction(
        intermediate,
        isolate_fastq_path,
        isolate_bam_path,
        logger,
        0.01,  # p_score_cutoff
        proc,
        results,
        run_subprocess,
        subtraction_indexes_path,
        subtractions,
        subtracted_bam_path,
        work_path,
    )

    assert results["subtracted_count"] == 0 if no_subtractions else 4

    def parse_headers(path: Path) -> dict:
        headers = {}

        with pysam.AlignmentFile(path, "r") as alignment_file:
            # Parse headers
            header = alignment_file.header.to_dict()

            # Process HD (header) entries
            if "HD" in header:
                headers["HD"] = header["HD"]

            # Process SQ (sequence) entries
            if "SQ" in header:
                headers["SQ"] = header["SQ"]

            # Process PG (program) entries, filtering out variable command lines
            if "PG" in header:
                filtered_pg = []
                for pg_entry in header["PG"]:
                    # Create a copy without the variable CL (command line) field
                    filtered_entry = {k: v for k, v in pg_entry.items() if k != "CL"}
                    filtered_pg.append(filtered_entry)

                headers["PG"] = filtered_pg

        return headers

    def parse_alignments(path: Path) -> set[tuple]:
        """Parse alignment file (SAM or BAM) into a standardized dict structure."""
        with pysam.AlignmentFile(path, "r") as alignment_file:
            return {
                (
                    read.query_name,
                    read.flag,
                    read.reference_name,
                    read.reference_start,
                    read.mapping_quality,
                    str(read.cigarstring) if read.cigarstring else "*",
                    read.template_length,
                    read.query_sequence[:5],
                )
                for read in alignment_file
            }

    assert parse_alignments(subtracted_bam_path) == snapshot(name="alignments")
    assert parse_headers(subtracted_bam_path) == snapshot(name="headers")


async def test_pathoscope(
    analysis: WFAnalysis,
    example_path: Path,
    index: WFIndex,
    mocker,
    snapshot: SnapshotAssertion,
    subtracted_bam_path: Path,
    work_path: Path,
):
    shutil.copyfile(example_path / "to_isolates.bam", subtracted_bam_path)

    intermediate = SimpleNamespace(candidate_sequence_ids={"candidate"})
    p_score_cutoff = 0.01
    results = {}

    logger = get_logger("test")

    # Mock upload_result to capture the hits
    upload_result_mock = mocker.patch.object(analysis, "upload_result")

    await reassignment(
        analysis,
        index,
        intermediate,
        logger,
        p_score_cutoff,
        results,
        subtracted_bam_path,
        work_path,
    )

    # Verify upload_result was called
    upload_result_mock.assert_called_once()

    # Extract the hits from the call
    call_args = upload_result_mock.call_args[0][0]
    hits = call_args["hits"]

    # Verify each hit has the expected structure
    for hit in hits:
        assert "id" in hit
        assert "otu" in hit
        assert "id" in hit["otu"]
        assert "version" in hit["otu"]
        assert "align" in hit
        assert "coverage" in hit
        assert "depth" in hit
        assert "final" in hit
        assert "initial" in hit
        assert isinstance(hit["align"], list)
        assert isinstance(hit["coverage"], float)
        assert isinstance(hit["depth"], (int, float))

    # Snapshot test for the hits structure
    assert hits == snapshot
    for hit in hits:
        if hit["id"].startswith("13TF149_Reovirus"):
            assert hit["otu"] == {"id": "reo", "version": 5}
        elif hit["id"] in {"JQ080272", "KX109927"}:
            assert hit["otu"] == {"id": "baz", "version": 6}
        else:
            assert hit["otu"] == {"id": "foobar", "version": 10}

    report: dict[str, list] = {}

    with open(work_path / "report.tsv") as f:  # noqa: ASYNC230
        for line in f:
            if "Final Guess" in line:
                continue

            if "Total Number of Aligned Reads" in line:
                continue

            split = line.split("\t")

            report[split[0]] = [float(f"{float(n):.5g}") for n in split[1:]]

    assert report == snapshot


async def test_build_isolate_index_resolves_candidate_sequences_to_otus(
    collapsed_reference_path: Path,
    index: WFIndex,
    isolate_index_path: Path,
    mocker,
    tmp_path: Path,
):
    await create_test_index(collapsed_reference_path)
    isolate_fasta_path = tmp_path / "isolates.fa"
    run_subprocess = mocker.AsyncMock()
    intermediate = SimpleNamespace(
        candidate_sequence_ids={"default-a", "non-default", "non-default-only"},
    )

    await build_isolate_index(
        collapsed_reference_path,
        index,
        intermediate,
        isolate_fasta_path,
        isolate_index_path,
        run_subprocess,
        2,
    )

    assert isolate_fasta_path.read_text().splitlines() == [
        ">default-a",
        "ACGT",
        ">default-b",
        "TTA",
        ">non-default",
        "GGGG",
        ">non-default-only",
        "CCCC",
    ]
    run_subprocess.assert_awaited_once_with(
        [
            "bowtie2-build",
            "--threads",
            "2",
            str(isolate_fasta_path),
            str(isolate_index_path),
        ],
    )


async def test_build_isolate_index_no_candidates(
    collapsed_reference_path: Path,
    index: WFIndex,
    mocker,
    work_path: Path,
):
    """When no candidate OTUs are found the isolate FASTA is empty.

    ``bowtie2-build`` exits 1 on an empty FASTA, so the index build must be skipped
    entirely (VIR-2569) rather than invoking the subprocess.
    """
    isolate_path = work_path / "isolates"
    isolate_path.mkdir()

    run_subprocess = mocker.AsyncMock()
    intermediate = SimpleNamespace(candidate_sequence_ids=set())

    await build_isolate_index(
        collapsed_reference_path,
        index,
        intermediate,
        isolate_path / "isolate_index.fa",
        isolate_path / "isolates",
        run_subprocess,
        2,
    )

    run_subprocess.assert_not_awaited()


async def test_no_candidates_uploads_empty_result(
    analysis: WFAnalysis,
    index: WFIndex,
    mocker,
    sample: WFSample,
    work_path: Path,
):
    """With no candidate OTUs the downstream steps short-circuit.

    ``map_isolates``, ``eliminate_subtraction`` and ``reassignment`` must not run any
    subprocess and the analysis is finished with an empty result (VIR-2569).
    """
    run_subprocess = mocker.AsyncMock()
    intermediate = SimpleNamespace(candidate_sequence_ids=set())
    logger = get_logger("test")
    results = {}

    await map_isolates(
        intermediate,
        work_path / "isolate_mapped.fq",
        work_path / "isolates",
        work_path / "to_isolates.bam",
        2,
        run_subprocess,
        sample,
    )

    await eliminate_subtraction(
        intermediate,
        work_path / "isolate_mapped.fq",
        work_path / "to_isolates.bam",
        logger,
        0.01,
        2,
        results,
        run_subprocess,
        work_path / "subtraction_indexes",
        [],
        work_path / "subtracted.bam",
        work_path,
    )

    assert results["subtracted_count"] == 0

    await reassignment(
        analysis,
        index,
        intermediate,
        logger,
        0.01,
        results,
        work_path / "subtracted.bam",
        work_path,
    )

    run_subprocess.assert_not_awaited()
    analysis.upload_result.assert_called_once_with(
        {"subtracted_count": 0, "read_count": 0, "hits": []},
    )
