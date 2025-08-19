import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
from structlog import get_logger
from syrupy import SnapshotAssertion
from virtool_workflow import RunSubprocess
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.indexes import WFIndex
from virtool_workflow.data.samples import WFSample
from virtool_workflow.data.subtractions import WFSubtraction
from virtool_workflow.pytest_plugin import Data

from workflow import (
    eliminate_subtraction,
    map_default_isolates,
    map_isolates,
    reassignment,
)


@pytest.fixture()
def work_path(tmpdir):
    path = Path(tmpdir) / "work"
    path.mkdir(parents=True)

    return path


@pytest.fixture()
def analysis(data: Data, mocker):
    analysis_ = mocker.Mock(WFAnalysis)

    analysis_.id = data.analysis.id
    analysis_.workflow = "pathoscope_bowtie"
    analysis_.ready = False
    analysis_.sample = data.analysis.sample

    return analysis_


@pytest.fixture()
def index(data: Data, example_path: Path, work_path: Path):
    data.index.manifest = {"foobar": 10, "reo": 5, "baz": 6}

    index_path = work_path / "indexes" / data.index.id

    shutil.copytree(example_path / "index", index_path)

    return WFIndex(
        id=data.index.id,
        path=index_path,
        manifest=data.index.manifest,
        reference=data.index.reference,
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
def sample(data: Data, example_path: Path, work_path: Path):
    data.sample.library_type = "normal"

    path = work_path / "samples" / data.sample.id
    path.mkdir(parents=True)

    shutil.copyfile(example_path / "sample" / "reads_1.fq.gz", path / "reads_1.fq.gz")

    return WFSample(
        id=data.sample.id,
        library_type=data.sample.library_type,
        name=data.sample.name,
        paired=False,
        quality=data.sample.quality,
        read_paths=(path / "reads_1.fq.gz",),
    )


@pytest.fixture()
def subtractions(data: Data, example_path: Path, work_path: Path):
    subtraction_path = work_path / "subtractions" / "subtraction"
    subtraction_path.parent.mkdir(parents=True)

    shutil.copytree(example_path / "subtraction", subtraction_path)

    return [
        WFSubtraction(
            id=data.subtraction.id,
            files=[],
            gc=data.subtraction.gc,
            name=data.subtraction.name,
            nickname=data.subtraction.nickname,
            path=subtraction_path,
        ),
    ]


async def test_map_default_isolates(
    index: WFIndex,
    run_subprocess,
    sample: WFSample,
    work_path: Path,
    snapshot: SnapshotAssertion,
):
    intermediate = SimpleNamespace(to_otus=set())

    logger = get_logger("test")

    await map_default_isolates(
        intermediate,
        logger,
        index,
        2,
        0.01,
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
    import pysam

    for path in (example_path / "index").iterdir():
        if "reference" in path.name:
            shutil.copyfile(
                path,
                work_path / path.name.replace("reference", "isolates"),
            )

    intermediate = SimpleNamespace(isolate_high_scores={})

    isolate_fastq_path = work_path / "mapped.fq"
    isolate_index_path = work_path / "isolates"
    isolate_bam_path = work_path / "to_isolates.bam"

    proc = 1
    p_score = 0.01

    await map_isolates(
        intermediate,
        isolate_fastq_path,
        isolate_index_path,
        isolate_bam_path,
        proc,
        p_score,
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

    # Note: isolate_high_scores is now populated during eliminate_subtraction step
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
    subtracted_path = work_path / "subtracted.sam"

    logger = get_logger("test")

    # Convert example SAM to BAM for this test
    await run_subprocess(
        [
            "samtools",
            "view",
            "-bS",
            "-o",
            str(isolate_bam_path),
            str(example_path / "to_isolates.sam"),
        ]
    )
    shutil.copyfile(example_path / "to_isolates.fq", isolate_fastq_path)

    proc = 2

    results = {}

    if no_subtractions:
        subtractions = []

    intermediate = SimpleNamespace()

    await eliminate_subtraction(
        intermediate,
        isolate_fastq_path,
        isolate_bam_path,
        logger,
        0.01,  # p_score_cutoff
        proc,
        results,
        run_subprocess,
        subtractions,
        subtracted_path,
        work_path,
    )

    assert results["subtracted_count"] == 0 if no_subtractions else 4

    assert not (work_path / "to_subtraction.sam").is_file()
    assert (work_path / "subtracted.sam").is_file()

    lines: dict[str, list] = {}

    with open(work_path / "subtracted.sam") as f:
        for line in f:
            split = line.split("\t")
            # Filter out @PG lines that contain variable paths
            if split[0].startswith("@PG") and len(split) > 1:
                filtered_parts = []
                for part in split[1:]:
                    # Skip the CL (command line) field that contains variable paths
                    if part.startswith("CL:"):
                        continue
                    filtered_parts.append(part)
                lines[split[0]] = filtered_parts
            else:
                lines[split[0]] = split[1:]

    assert lines == snapshot


async def test_pathoscope(
    analysis: WFAnalysis,
    example_path: Path,
    index: WFIndex,
    mocker,
    ref_lengths,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    subtracted_sam_path = work_path / "subtracted.sam"
    shutil.copyfile(example_path / "to_isolates.sam", subtracted_sam_path)

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
        subtracted_sam_path,
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
