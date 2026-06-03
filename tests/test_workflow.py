import gzip
import json
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
from virtool.workflow.data.cache import CacheHit, CacheMiss
from virtool.workflow.data.indexes import WFIndex
from virtool.workflow.data.samples import WFSample
from virtool.workflow.data.subtractions import WFSubtraction
from virtool.workflow.pytest_plugin import WorkflowData

from workflow import (
    create_reference_index,
    create_subtraction_index,
    eliminate_subtraction,
    get_mapping_index_cache_params,
    map_default_isolates,
    map_isolates,
    reassignment,
)
from workflow_pathoscope.utils import write_default_isolate_fasta


BOWTIE2_INDEX_SUFFIXES = (
    "1.bt2",
    "2.bt2",
    "3.bt2",
    "4.bt2",
    "rev.1.bt2",
    "rev.2.bt2",
)


@pytest.fixture()
def work_path(tmpdir):
    path = Path(tmpdir) / "work"
    path.mkdir(parents=True)

    return path


@pytest.fixture()
def reference_index_path(work_path: Path) -> Path:
    return work_path / "reference_index" / "reference"


@pytest.fixture()
def subtraction_index_path(
    subtractions: list[WFSubtraction],
    work_path: Path,
) -> Path:
    return work_path / "subtraction_indexes" / subtractions[0].id / "subtraction"


class FakeWorkflowCache:
    def __init__(self, hit_source: Path | None = None):
        self.hit_source = hit_source
        self.gets = []
        self.puts = []

    async def get(self, key: str, target: Path):
        self.gets.append((key, target))

        if self.hit_source is None:
            return CacheMiss(key)

        target.mkdir(parents=True, exist_ok=True)

        for path in self.hit_source.iterdir():
            shutil.copyfile(path, target / path.name)

        return CacheHit(key, target)

    async def put(self, key: str, source: Path, params: dict | None = None):
        self.puts.append((key, source, params))
        return True


class FakeRunSubprocess:
    def __init__(self):
        self.commands = []

    async def __call__(
        self,
        command: list[str],
        cwd: str | Path | None = None,
        env: dict | None = None,
        stderr_handler=None,
        stdout_handler=None,
    ):
        self.commands.append(command)

        if command == ["bowtie2-build", "--version"]:
            await stdout_handler(b"/usr/bin/bowtie2-build-s version 2.5.4\n")
            return SimpleNamespace(returncode=0)

        if command[0] == "bowtie2-build":
            prefix = Path(command[-1])
            prefix.parent.mkdir(parents=True, exist_ok=True)

            for suffix in BOWTIE2_INDEX_SUFFIXES:
                (prefix.parent / f"{prefix.name}.{suffix}").write_bytes(
                    f"{prefix.name}.{suffix}".encode(),
                )

            return SimpleNamespace(returncode=0)

        raise AssertionError(f"Unexpected subprocess command: {command}")


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
                                "default": True,
                                "sequences": [
                                    {"_id": "default-a", "sequence": "ACGT"},
                                    {"_id": "default-b", "sequence": "TTA"},
                                ],
                            },
                            {
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
    index: WFIndex,
    reference_index_path: Path,
    tmp_path: Path,
):
    source = tmp_path / "cached-reference"
    write_bowtie2_bundle(
        source,
        "reference",
        b"cached-reference",
        {"reference.fa": b"cached-fasta", "cache-manifest.json": b"cached-manifest"},
    )
    cache = FakeWorkflowCache(source)
    run_subprocess = FakeRunSubprocess()
    logger = get_logger("test")

    result_index_path = await create_reference_index(
        cache,
        index,
        logger,
        4,
        run_subprocess,
        reference_index_path,
    )

    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        FakeRunSubprocess(),
        {
            "source": "index.json_path",
            "selection": "default_isolates",
        },
    )

    assert cache.gets == [
        (
            derive_key(params),
            reference_index_path.parent,
        ),
    ]
    assert cache.puts == []
    assert result_index_path == reference_index_path
    assert run_subprocess.commands == [["bowtie2-build", "--version"]]
    assert read_directory_bytes(reference_index_path.parent) == read_directory_bytes(
        source
    )


async def test_create_reference_index_miss(
    index: WFIndex,
    reference_index_path: Path,
):
    write_reference_json(index.json_path)
    cache = FakeWorkflowCache()
    run_subprocess = FakeRunSubprocess()
    logger = get_logger("test")

    result_index_path = await create_reference_index(
        cache,
        index,
        logger,
        4,
        run_subprocess,
        reference_index_path,
    )

    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        FakeRunSubprocess(),
        {
            "source": "index.json_path",
            "selection": "default_isolates",
        },
    )
    key = derive_key(params)
    assert cache.gets == [(key, reference_index_path.parent)]
    assert cache.puts == [(key, reference_index_path.parent, params)]
    assert result_index_path == reference_index_path
    assert params == {
        "index_kind": "reference_mapping_index",
        "workflow": "pathoscope",
        "parent_id": index.id,
        "source": "index.json_path",
        "selection": "default_isolates",
        "tool_name": "bowtie2-build",
        "tool_version": "2.5.4",
    }
    assert (reference_index_path.parent / "reference.fa").read_text() == (
        ">default-a\nACGT\n>default-b\nTTA\n"
    )
    assert run_subprocess.commands == [
        ["bowtie2-build", "--version"],
        [
            "bowtie2-build",
            "--threads",
            "4",
            str(reference_index_path.parent / "reference.fa"),
            str(reference_index_path),
        ],
    ]
    assert read_directory_bytes(reference_index_path.parent) == {
        "reference.fa": b">default-a\nACGT\n>default-b\nTTA\n",
        **bowtie2_bundle_bytes("reference"),
    }


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
    subtractions: list[WFSubtraction],
    subtraction_index_path: Path,
    tmp_path: Path,
    work_path: Path,
):
    source = tmp_path / "cached-subtraction"
    write_bowtie2_bundle(
        source,
        "subtraction",
        b"cached-subtraction",
        {"cache-manifest.json": b"cached-manifest"},
    )
    cache = FakeWorkflowCache(source)
    run_subprocess = FakeRunSubprocess()
    logger = get_logger("test")
    subtraction = subtractions[0]
    subtraction_index_paths = {}

    result_index_paths = await create_subtraction_index(
        cache,
        logger,
        4,
        run_subprocess,
        subtractions,
        subtraction_index_paths,
        work_path,
    )

    params = await get_mapping_index_cache_params(
        "subtraction_mapping_index",
        subtraction.id,
        FakeRunSubprocess(),
    )

    assert cache.gets == [
        (
            derive_key(params),
            subtraction_index_path.parent,
        ),
    ]
    assert cache.puts == []
    assert result_index_paths == {subtraction.id: subtraction_index_path}
    assert subtraction_index_paths == {subtraction.id: subtraction_index_path}
    assert run_subprocess.commands == [["bowtie2-build", "--version"]]
    assert read_directory_bytes(subtraction_index_path.parent) == read_directory_bytes(
        source
    )


async def test_create_subtraction_index_miss(
    subtractions: list[WFSubtraction],
    subtraction_index_path: Path,
    work_path: Path,
):
    cache = FakeWorkflowCache()
    run_subprocess = FakeRunSubprocess()
    logger = get_logger("test")
    subtraction = subtractions[0]
    subtraction_index_paths = {}

    result_index_paths = await create_subtraction_index(
        cache,
        logger,
        4,
        run_subprocess,
        subtractions,
        subtraction_index_paths,
        work_path,
    )

    params = await get_mapping_index_cache_params(
        "subtraction_mapping_index",
        subtraction.id,
        FakeRunSubprocess(),
    )
    key = derive_key(params)
    assert cache.gets == [(key, subtraction_index_path.parent)]
    assert cache.puts == [(key, subtraction_index_path.parent, params)]
    assert result_index_paths == {subtraction.id: subtraction_index_path}
    assert subtraction_index_paths == {subtraction.id: subtraction_index_path}
    assert params == {
        "index_kind": "subtraction_mapping_index",
        "workflow": "pathoscope",
        "parent_id": subtraction.id,
        "tool_name": "bowtie2-build",
        "tool_version": "2.5.4",
    }
    assert run_subprocess.commands == [
        ["bowtie2-build", "--version"],
        [
            "bowtie2-build",
            "--threads",
            "4",
            str(subtraction.fasta_path),
            str(subtraction_index_path),
        ],
    ]
    assert read_directory_bytes(subtraction_index_path.parent) == bowtie2_bundle_bytes(
        "subtraction"
    )


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
    sample: WFSample,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    for path in (example_path / "index").iterdir():
        if "reference" in path.name:
            shutil.copyfile(
                path,
                work_path / path.name.replace("reference", "isolates"),
            )

    isolate_fastq_path = work_path / "mapped.fq"
    isolate_index_path = work_path / "isolates"
    isolate_bam_path = work_path / "to_isolates.bam"

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
    no_subtractions: bool,
    subtractions: list[WFSubtraction],
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    isolate_fastq_path = work_path / "to_isolates.fq"
    isolate_bam_path = work_path / "to_isolates.bam"
    subtracted_path = work_path / "subtracted.bam"

    shutil.copyfile(example_path / "to_isolates.bam", isolate_bam_path)
    shutil.copyfile(example_path / "to_isolates.fq", isolate_fastq_path)

    logger = get_logger("test")
    proc = 2
    results = {}

    if no_subtractions:
        subtractions = []

    intermediate = SimpleNamespace()
    subtraction_index_paths = {}

    if subtractions:
        await create_subtraction_index(
            FakeWorkflowCache(subtractions[0].path),
            logger,
            proc,
            FakeRunSubprocess(),
            subtractions,
            subtraction_index_paths,
            work_path,
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
        subtraction_index_paths,
        subtractions,
        subtracted_path,
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

    assert parse_alignments(work_path / "subtracted.bam") == snapshot(name="alignments")
    assert parse_headers(work_path / "subtracted.bam") == snapshot(name="headers")


async def test_pathoscope(
    analysis: WFAnalysis,
    example_path: Path,
    index: WFIndex,
    mocker,
    ref_lengths,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    subtracted_bam_path = work_path / "to_isolates.bam"

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
