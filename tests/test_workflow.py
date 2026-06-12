import asyncio
import gzip
import json
import re
import shutil

import pysam
from pathlib import Path
from types import SimpleNamespace

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
    collapse_reference,
    create_reference_index,
    create_subtraction_index,
    eliminate_subtraction,
    get_subtraction_index_path,
    map_default_isolates,
    map_isolates,
    reassignment,
)
from workflow_pathoscope.utils import (
    CD_HIT_EST_IDENTITY,
    collapse_reference_json,
    get_mapping_index_cache_params,
    get_reference_collapse_cache_params,
    write_default_isolate_fasta,
    write_isolate_fasta,
)


BOWTIE2_INDEX_SUFFIXES = (
    "1.bt2",
    "2.bt2",
    "3.bt2",
    "4.bt2",
    "rev.1.bt2",
    "rev.2.bt2",
)
REDUNDANT_REFERENCE_JSON_PATH = (
    Path(__file__).parent / "assets" / "redundant_reference.json"
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


def write_reference_json(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "wt") as f:
        json.dump(
            {
                "otus": [
                    {
                        "_id": "default-otu",
                        "isolates": [
                            {
                                "id": "default",
                                "default": True,
                                "sequences": [
                                    {"_id": "default-a", "sequence": "ACGT"},
                                    {"_id": "default-b", "sequence": "TTA"},
                                ],
                            },
                            {
                                "id": "non-default",
                                "default": False,
                                "sequences": [
                                    {"_id": "non-default", "sequence": "GGGG"},
                                ],
                            },
                        ],
                    },
                    {
                        "_id": "non-default-otu",
                        "isolates": [
                            {
                                "id": "non-default-only",
                                "default": False,
                                "sequences": [
                                    {"_id": "non-default-only", "sequence": "CCCC"},
                                ],
                            },
                        ],
                    },
                ],
            },
            f,
        )


def write_redundant_reference_json(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(REDUNDANT_REFERENCE_JSON_PATH, path)


@pytest.fixture()
def analysis(workflow_data: WorkflowData, mocker):
    analysis_ = mocker.Mock(WFAnalysis)

    analysis_.id = workflow_data.analysis.id
    analysis_.workflow = "pathoscope_bowtie"
    analysis_.ready = False
    analysis_.sample = workflow_data.analysis.sample

    return analysis_


@pytest.fixture()
def index(workflow_data: WorkflowData, example_path: Path, work_path: Path):
    workflow_data.index.manifest = {"foobar": 10, "reo": 5, "baz": 6}

    index_path = work_path / "indexes" / workflow_data.index.id

    shutil.copytree(example_path / "index", index_path)

    return WFIndex(
        id=workflow_data.index.id,
        path=index_path,
        manifest=workflow_data.index.manifest,
        reference=workflow_data.index.reference,
        sequence_lengths={},
        sequence_otu_map={
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
        },
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
    write_reference_json(source / collapsed_reference_path.name)
    (source / "collapse-manifest.json").write_text(
        json.dumps(
            {
                "isolate_count_before": 4,
                "isolate_count_after": 3,
                "isolate_count_removed": 1,
            }
        )
    )
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
    index: WFIndex,
    run_subprocess: RunSubprocess,
    workflow_cache: WorkflowCache,
):
    write_redundant_reference_json(index.json_path)
    logger = get_logger("test")

    result_path = await collapse_reference(
        workflow_cache,
        collapsed_reference_path,
        index,
        logger,
        4,
        run_subprocess,
    )

    params = await get_reference_collapse_cache_params(index.id, run_subprocess)

    assert result_path == collapsed_reference_path
    assert_cache_params(
        params,
        {
            "index_kind": "collapsed_reference",
            "workflow": "pathoscope",
            "workflow_version": "UNKNOWN",
            "parent_id": index.id,
            "source": "index_json",
            "tool_name": "cd-hit-est",
            "identity": "0.99",
        },
    )

    with open(collapsed_reference_path) as handle:
        collapsed_reference = json.load(handle)

    assert [
        isolate["id"] for isolate in collapsed_reference["otus"][0]["isolates"]
    ] == [
        "default",
        "representative-1",
        "representative-2",
        "unique-combo",
    ]
    assert [
        (isolate["sequences"][0]["_id"], isolate["sequences"][1]["_id"])
        for isolate in collapsed_reference["otus"][0]["isolates"]
    ] == [
        ("default-a", "default-b"),
        ("representative-1-a", "representative-1-b"),
        ("representative-2-a", "representative-2-b"),
        ("unique-combo-a", "unique-combo-b"),
    ]
    assert json.loads(
        (collapsed_reference_path.parent / "collapse-manifest.json").read_text()
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 4,
        "isolate_count_removed": 1,
    }


async def test_collapse_reference_json_outputs_collapsed_reference_fasta(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    source_path = tmp_path / "reference.json"
    collapsed_path = tmp_path / "collapsed" / "reference.json"
    default_fasta_path = tmp_path / "default.fa"
    isolate_fasta_path = tmp_path / "isolates.fa"

    write_redundant_reference_json(source_path)

    assert await collapse_reference_json(
        source_path,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 4,
        "isolate_count_removed": 1,
    }

    assert write_default_isolate_fasta(collapsed_path, default_fasta_path) == {
        "default-a": 20,
        "default-b": 20,
    }
    assert write_isolate_fasta(
        {"collapse-otu"},
        collapsed_path,
        isolate_fasta_path,
    ) == {
        "default-a": 20,
        "default-b": 20,
        "representative-1-a": 20,
        "representative-1-b": 20,
        "representative-2-a": 20,
        "representative-2-b": 20,
        "unique-combo-a": 20,
        "unique-combo-b": 20,
    }
    assert "duplicate-a" not in isolate_fasta_path.read_text()
    assert "duplicate-b" not in isolate_fasta_path.read_text()


async def test_collapse_reference_json_allows_isolates_missing_schema_segments(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    source_path = tmp_path / "reference.json"
    collapsed_path = tmp_path / "collapsed" / "reference.json"

    write_redundant_reference_json(source_path)

    with open(source_path) as handle:
        reference_data = json.load(handle)

    reference_data["otus"][0]["isolates"][0]["sequences"] = [
        reference_data["otus"][0]["isolates"][0]["sequences"][1],
    ]

    with open(source_path, "w") as handle:
        json.dump(reference_data, handle)

    assert await collapse_reference_json(
        source_path,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 4,
        "isolate_count_removed": 1,
    }

    with open(collapsed_path) as handle:
        collapsed_reference = json.load(handle)

    assert collapsed_reference["otus"][0]["isolates"][0]["sequences"] == [
        {
            "_id": "default-b",
            "segment": "b",
            "sequence": "TGCATGCATGCATGCATGCA",
        },
    ]


async def test_collapse_reference_json_allows_unsegmented_isolates(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    source_path = tmp_path / "reference.json"
    collapsed_path = tmp_path / "collapsed" / "reference.json"

    write_redundant_reference_json(source_path)

    with open(source_path) as handle:
        reference_data = json.load(handle)

    reference_data["otus"][0]["schema"] = []

    for isolate in reference_data["otus"][0]["isolates"]:
        for sequence in isolate["sequences"]:
            sequence["segment"] = ""

    with open(source_path, "w") as handle:
        json.dump(reference_data, handle)

    assert await collapse_reference_json(
        source_path,
        collapsed_path,
        2,
        run_subprocess,
    ) == {
        "isolate_count_before": 5,
        "isolate_count_after": 4,
        "isolate_count_removed": 1,
    }

    with open(collapsed_path) as handle:
        collapsed_reference = json.load(handle)

    assert collapsed_reference["otus"][0]["schema"] == []
    assert {
        sequence["segment"]
        for isolate in collapsed_reference["otus"][0]["isolates"]
        for sequence in isolate["sequences"]
    } == {""}


async def test_collapse_reference_json_rejects_isolate_segments_that_do_not_match_schema(
    run_subprocess: RunSubprocess,
    tmp_path: Path,
):
    source_path = tmp_path / "reference.json"
    collapsed_path = tmp_path / "collapsed" / "reference.json"

    write_redundant_reference_json(source_path)

    with open(source_path) as handle:
        reference_data = json.load(handle)

    reference_data["otus"][0]["isolates"][0]["sequences"][0]["segment"] = "c"

    with open(source_path, "w") as handle:
        json.dump(reference_data, handle)

    with pytest.raises(
        ValueError,
        match=re.escape(
            "Isolate default sequence segments ['c'] are not defined in OTU collapse-otu "
            "schema segments ['a', 'b']"
        ),
    ):
        await collapse_reference_json(
            source_path,
            collapsed_path,
            2,
            run_subprocess,
        )


@pytest.fixture()
def sample(workflow_data: WorkflowData, example_path: Path, work_path: Path):
    workflow_data.sample.library_type = "normal"

    path = work_path / "samples" / workflow_data.sample.id
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
    write_reference_json(collapsed_reference_path)
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
    write_reference_json(collapsed_reference_path)
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
    write_reference_json(collapsed_reference_path)
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
    write_reference_json(collapsed_reference_path)
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


def test_write_default_isolate_fasta(tmp_path: Path):
    json_path = tmp_path / "reference.json.gz"
    target_path = tmp_path / "reference.fa"

    write_reference_json(json_path)

    assert write_default_isolate_fasta(json_path, target_path) == {
        "default-a": 4,
        "default-b": 3,
    }
    assert target_path.read_text() == ">default-a\nACGT\n>default-b\nTTA\n"


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

    intermediate = SimpleNamespace(to_otus=set())

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

    assert sorted(intermediate.to_otus) == snapshot


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

    await map_isolates(
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

    intermediate = SimpleNamespace()

    if subtractions:
        cached_subtraction_path = tmp_path / subtractions[0].id
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
    ref_lengths,
    snapshot: SnapshotAssertion,
    subtracted_bam_path: Path,
    work_path: Path,
):
    shutil.copyfile(example_path / "to_isolates.bam", subtracted_bam_path)

    intermediate = SimpleNamespace(lengths=ref_lengths)
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

    report: dict[str, list] = {}

    with open(work_path / "report.tsv") as f:
        for line in f:
            if "Final Guess" in line:
                continue

            if "Total Number of Aligned Reads" in line:
                continue

            split = line.split("\t")

            report[split[0]] = [float(f"{float(n):.5g}") for n in split[1:]]

    assert report == snapshot
