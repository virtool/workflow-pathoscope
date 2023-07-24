import json
import shutil
from pathlib import Path
from shutil import copytree
from types import SimpleNamespace

import pytest
from aiohttp.test_utils import make_mocked_coro
from pydantic_factories import ModelFactory
from virtool_core.models.index import Index
from virtool_core.models.subtraction import (
    Subtraction,
)
from virtool_workflow.data_model.analysis import WFAnalysis
from virtool_workflow.data_model.indexes import WFIndex
from virtool_workflow.data_model.samples import Sample, WFSample
from virtool_workflow.data_model.subtractions import WFSubtraction
from virtool_workflow.runtime.run_subprocess import run_subprocess as wf_run_subprocess

from workflow import (
    eliminate_subtraction,
    map_default_isolates,
    map_isolates,
    reassignment,
)

TEST_DATA_PATH = Path(__file__).parent / "test_files"
FASTQ_PATH = TEST_DATA_PATH / "test.fq"
INDEX_PATH = TEST_DATA_PATH / "index"
SUBTRACTION_PATH = TEST_DATA_PATH / "subtraction"
SAM_PATH = TEST_DATA_PATH / "test_al.sam"


@pytest.fixture
def work_path(tmpdir):
    return Path(tmpdir)


@pytest.fixture
def index(work_path: Path):
    index_path = work_path / "indexes/index3"
    shutil.copytree(INDEX_PATH, index_path)

    class IndexFactory(ModelFactory):
        __model__ = Index

    manifest = {"foobar": 10, "reo": 5, "baz": 6}

    index = WFIndex(
        IndexFactory.build(manifest=manifest),
        index_path,
        make_mocked_coro(),
        make_mocked_coro(),
        run_subprocess,
    )

    index._sequence_otu_map = {
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

    return index


@pytest.fixture
def sample(work_path: Path):
    shutil.copyfile(FASTQ_PATH, work_path / "reads_1.fq.gz")

    class SampleFactory(ModelFactory):
        __model__ = Sample

    _sample = WFSample.parse_obj(SampleFactory.build())

    _sample.read_paths = (work_path / "reads_1.fq.gz",)

    return _sample


@pytest.fixture
def read_file_names(sample: WFSample):
    return ",".join(str(p) for p in sample.read_paths)


@pytest.fixture
def subtractions(work_path: Path):
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    subtraction_path = work_path / "subtractions" / "subtraction"

    copytree(SUBTRACTION_PATH, subtraction_path)

    class SubtractionFactory(ModelFactory):
        __model__ = Subtraction

    _subtraction1 = SubtractionFactory.build()
    _subtraction2 = SubtractionFactory.build()
    _subtraction3 = SubtractionFactory.build()

    return [
        WFSubtraction(
            **{**_subtraction1.dict(), "ready": True},
            path=subtraction_path,
        ),
        WFSubtraction(**{**_subtraction2.dict(), "ready": True}, path=subtraction_path),
        WFSubtraction(**{**_subtraction3.dict(), "ready": True}, path=subtraction_path),
    ]


@pytest.fixture
async def run_subprocess():
    return wf_run_subprocess()


@pytest.fixture
def ref_lengths():
    with (TEST_DATA_PATH / "ref_lengths.json").open("r") as f:
        return json.load(f)


async def test_map_default_isolates(
    read_file_names, index: WFIndex, run_subprocess, snapshot
):
    intermediate = SimpleNamespace(to_otus=set())

    await map_default_isolates(
        intermediate, read_file_names, index, 2, 0.01, run_subprocess
    )

    assert sorted(intermediate.to_otus) == snapshot


async def test_map_isolates(
    index,
    read_file_names,
    work_path,
    run_subprocess,
    snapshot,
):
    for path in INDEX_PATH.iterdir():
        if "reference" in path.name:
            shutil.copyfile(
                path, work_path / path.name.replace("reference", "isolates")
            )

    intermediate = SimpleNamespace(isolate_high_scores={})
    isolate_fastq_path = work_path / "mapped.fq"
    isolate_sam_path = work_path / "to_isolates.sam"

    await map_isolates(
        read_file_names,
        intermediate,
        isolate_fastq_path,
        work_path / "isolates",
        isolate_sam_path,
        run_subprocess,
        1,
        0.01,
    )

    with isolate_sam_path.open("r") as f:
        assert sorted(line.rstrip() for line in f) == snapshot

    assert intermediate.isolate_high_scores == snapshot


@pytest.mark.parametrize("none", [True, False])
@pytest.mark.datafiles(SAM_PATH, FASTQ_PATH)
async def test_eliminate_subtraction(
    none, datafiles, subtractions, work_path, run_subprocess
):
    isolate_fastq_path = work_path / "test.fq"
    isolate_sam_path = work_path / "test_al.sam"

    subtracted_path = work_path / "subtracted.sam"

    results = {}

    if none:
        subtractions = []

    await eliminate_subtraction(
        isolate_fastq_path,
        isolate_sam_path,
        2,
        results,
        run_subprocess,
        subtractions,
        subtracted_path,
        work_path,
    )

    if none:
        assert results["subtracted_count"] == 0
    else:
        assert results["subtracted_count"] == 4

    assert not (work_path / "to_subtraction.sam").is_file()
    assert (work_path / "subtracted.sam").is_file()


async def test_pathoscope(
    mocker,
    data_regression,
    file_regression,
    index,
    ref_lengths,
    snapshot,
    work_path: Path,
):
    subtracted_path = work_path / "subtracted.sam"
    shutil.copyfile(SAM_PATH, subtracted_path)

    analysis = mocker.Mock(spec=WFAnalysis)
    intermediate = SimpleNamespace(lengths=ref_lengths)
    results = {}

    await reassignment(
        analysis,
        index,
        intermediate,
        0.01,
        results,
        subtracted_path,
        work_path,
    )

    # check run function produces desired output
    with open(
        Path(__file__).parent / "test_workflow/test_pathoscope.tsv", "r"
    ) as testFile:
        for line in testFile.readlines()[2:]:
            newLine = line.split("\t")
            records = list(results["hits"])
            for record in records:
                if newLine[0] == record["id"]:
                    assert int(float(newLine[1]) - record["final"]["pi"]) == 0
                    assert int(float(newLine[2]) - record["final"]["best"]) == 0
                    assert int(float(newLine[4]) - record["final"]["high"]) == 0
                    assert int(float(newLine[5]) - record["final"]["low"]) == 0
                    assert int(float(newLine[3]) - record["final"]["reads"]) == 0
