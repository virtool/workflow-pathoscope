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
        run_subprocess,
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

    intermediate = SimpleNamespace(isolate_high_scores={})

    isolate_fastq_path = work_path / "mapped.fq"
    isolate_index_path = work_path / "isolates"
    isolate_sam_path = work_path / "to_isolates.sam"

    proc = 1
    p_score = 0.01

    await map_isolates(
        intermediate,
        isolate_fastq_path,
        isolate_index_path,
        isolate_sam_path,
        proc,
        p_score,
        run_subprocess,
        sample,
    )

    with isolate_sam_path.open("r") as f:
        assert sorted(line.rstrip() for line in f) == snapshot

    assert intermediate.isolate_high_scores == snapshot


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
    isolate_sam_path = work_path / "to_isolates.sam"
    subtracted_path = work_path / "subtracted.sam"

    logger = get_logger("test")

    shutil.copyfile(example_path / "to_isolates.sam", isolate_sam_path)
    shutil.copyfile(example_path / "to_isolates.fq", isolate_fastq_path)

    proc = 2

    results = {}

    if no_subtractions:
        subtractions = []

    await eliminate_subtraction(
        isolate_fastq_path,
        isolate_sam_path,
        logger,
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
            lines[split[0]] = split[1:]

    assert lines == snapshot


async def test_pathoscope(
    analysis: WFAnalysis,
    example_path: Path,
    index: WFIndex,
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
