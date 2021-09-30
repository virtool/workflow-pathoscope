import json
import logging
import shutil
from pathlib import Path
from shutil import copytree
from types import SimpleNamespace

import pytest
import virtool_workflow.execution.run_subprocess
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.library_types import LibraryType
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.data_model import NucleotideComposition, Subtraction
from virtool_workflow.data_model.samples import Sample

from workflow import map_default_isolates, map_subtractions, map_isolates, subtract_mapping, reassignment

logger = logging.getLogger(__name__)

TEST_DATA_PATH = Path(__file__).parent / "test_files"
FASTQ_PATH = TEST_DATA_PATH / "test.fq"
INDEX_PATH = TEST_DATA_PATH / "index"
REF_LENGTHS_PATH = TEST_DATA_PATH / "ref_lengths.json"
SUBTRACTION_PATH = TEST_DATA_PATH / "subtraction"
TO_SUBTRACTION_PATH = TEST_DATA_PATH / "to_subtraction.json"
VTA_PATH = TEST_DATA_PATH / "test.vta"


@pytest.fixture
def work_path(tmpdir):
    return Path(tmpdir)


@pytest.fixture
def index(work_path: Path):
    index_path = work_path / "indexes/index3"
    shutil.copytree(INDEX_PATH, index_path)

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
        "NC_007448": "foobar"
    }

    return Index(
        id="index3",
        manifest={
            "foobar": 10,
            "reo": 5,
            "baz": 6,
        },
        _sequence_otu_map=sequence_otu_map,
        reference=None,
        path=index_path,
        ready=True,
        _run_in_executor=run_in_executor,
        _run_subprocess=run_subprocess,
    )


@pytest.fixture
def sample():
    return Sample(
        id="foobar",
        name="foobar",
        paired=False,
        library_type=LibraryType.other,
        quality={
            "count": 1337,
            "length": [78, 101]
        },
        locale=None,
        isolate=None,
        host=None,
    )


@pytest.fixture
def reads(work_path: Path):
    shutil.copyfile(FASTQ_PATH, work_path / "reads_1.fq.gz")

    return Reads(
        sample=sample,
        quality={},
        path=work_path
    )


@pytest.fixture
def subtraction(work_path: Path):
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    subtraction_path = work_path / "subtractions" / "subtraction"

    copytree(SUBTRACTION_PATH, subtraction_path)

    nucleotide_composition = NucleotideComposition(
        a=0.1,
        t=0.2,
        g=0.3,
        c=0.4,
    )

    return Subtraction(
        id="arabidopsis_thaliana",
        name="Arabidopsis thaliana",
        nickname="Thalecress",
        count=12,
        gc=nucleotide_composition,
        path=subtraction_path,
    )


@pytest.fixture
async def run_in_executor():
    async def _run_in_executor(func, *args):
        return func(*args)

    return _run_in_executor


@pytest.fixture
async def run_subprocess():
    return virtool_workflow.execution.run_subprocess.run_subprocess()


@pytest.fixture
def ref_lengths():
    with REF_LENGTHS_PATH.open("r") as f:
        return json.load(f)


async def test_map_default_isolates(data_regression, reads: Reads, index: Index, run_subprocess):
    intermediate = SimpleNamespace(to_otus=set())

    await map_default_isolates(
        intermediate,
        reads,
        index,
        2,
        0.01,
        run_subprocess,
        logger
    )

    assert sorted(intermediate.to_otus) == sorted([
        "NC_013110",
        "NC_017938",
        "NC_006057",
        "NC_007448",
        "JQ080272",
        "NC_001836",
        "NC_003347",
        "NC_016509",
        "NC_017939",
        "NC_006056",
        "NC_003623",
        "KX109927",
        "NC_016416",
        "NC_001948",
        "NC_021148",
        "NC_003615",
        "NC_004006"
    ])


async def test_map_isolates(
        data_regression,
        index,
        reads: Reads,
        work_path,
        run_subprocess,
):
    for path in INDEX_PATH.iterdir():
        if "reference" in path.name:
            logger.info(path)
            shutil.copyfile(
                path,
                work_path / path.name.replace("reference", "isolates")
            )

    isolate_fastq_path = work_path / "mapped.fq"
    isolate_vta_path = work_path / "isolates.vta"

    await map_isolates(
        reads,
        isolate_fastq_path,
        work_path / "isolates",
        isolate_vta_path,
        run_subprocess,
        1,
        0.01
    )

    with isolate_vta_path.open("r") as f:
        data_regression.check(sorted([line.rstrip() for line in f]))


async def test_map_subtraction(
        data_regression,
        subtraction,
        work_path,
        run_subprocess
):
    isolate_fastq_path = work_path / "mapped.fastq"
    shutil.copyfile(FASTQ_PATH, isolate_fastq_path)

    intermediate = SimpleNamespace(
        isolate_mapped_fastq_path=isolate_fastq_path,
    )

    await map_subtractions(
        intermediate,
        subtraction,
        run_subprocess,
        2
    )

    data_regression.check(sorted(intermediate.to_subtraction))

@pytest.mark.datafiles(VTA_PATH)
async def test_subtract_mapping(
        file_regression,
        datafiles,
        run_in_executor,
        work_path: Path
):
    with TO_SUBTRACTION_PATH.open("r") as f:
        intermediate = SimpleNamespace(
            to_subtraction=json.load(f)
        )

    results = dict()

    await subtract_mapping(
        intermediate,
        datafiles / "test.vta",
        results,
        run_in_executor,
        work_path
    )

    assert results["subtracted_count"] == 4

    with open(work_path / "test.vta", "r") as f:
        file_regression.check(f.read())


async def test_pathoscope(
        data_regression,
        index,
        ref_lengths,
        run_in_executor,
        work_path: Path
):
    isolate_vta_path = work_path / "to_isolates.vta"

    intermediate = SimpleNamespace(
        lengths=ref_lengths
    )

    results = dict()

    shutil.copyfile(VTA_PATH, isolate_vta_path)

    await reassignment(
        intermediate,
        results,
        run_in_executor,
        isolate_vta_path,
        index,
        0.01,
        work_path
    )

    with intermediate.reassigned_path.open("r") as f:
        reassigned_data = sorted([line.rstrip() for line in f])

    with intermediate.report_path.open("r") as f:
        report_data = sorted([line.rstrip() for line in f])

    data_regression.check([
        reassigned_data,
        report_data,
        results
    ])
