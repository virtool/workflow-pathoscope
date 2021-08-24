import json
import logging
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
from fixtures import FixtureScope
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.library_types import LibraryType
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.data_model.samples import Sample

import workflow

logger = logging.getLogger(__name__)


@pytest.fixture
def sequence_otu_map():
    return {
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


@pytest.fixture
async def scope(
    tmpdir,
    fastq_path,
    otu_resource,
    indexes_path,
    sequence_otu_map,
):
    """Fixture scope for pathoscope workflow tests"""
    async with FixtureScope() as scope_:
        # Create a mock sample fixture
        scope_["sample"] = Sample(
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

        # Create a mock reads fixture
        temp_read_path = Path(tmpdir) / "reads"
        temp_read_path.mkdir()

        shutil.copyfile(fastq_path, temp_read_path/"reads_1.fq.gz")

        scope_["reads"] = Reads(
            sample=scope_["sample"],
            quality={},
            path=temp_read_path
        )

        run_in_executor = await scope_.get_or_instantiate("run_in_executor")
        run_subprocess = await scope_.get_or_instantiate("run_subprocess")

        # Create a mock index fixture
        temp_index_path = Path(tmpdir) / "indexes"
        shutil.copytree(indexes_path, temp_index_path)

        scope_["index"] = Index(
            id="index3",
            manifest={
                "foobar": 10,
                "reo": 5,
                "baz": 6,
            },
            _sequence_otu_map=sequence_otu_map,
            reference=None,
            path=temp_index_path,
            ready=True,
            _run_in_executor=run_in_executor,
            _run_subprocess=run_subprocess,
        )

        scope_["logger"] = logger

        yield scope_


@pytest.fixture
def otu_resource(vta_path):
    map_dict = dict()
    otus = dict()

    with open(vta_path, "r") as handle:
        for line in handle:
            ref_id = line.split(",")[1]

            otu_id = "otu_{}".format(ref_id)

            map_dict[ref_id] = otu_id

            otus[otu_id] = {
                "otu": otu_id,
                "version": 2
            }

    return map_dict, otus


async def test_map_default_isolates(tmpdir, fastq_path, scope: FixtureScope):
    scope["intermediate"] = SimpleNamespace(to_otus=set())

    map_default_isolates = await scope.bind(workflow.map_default_isolates)

    await map_default_isolates()

    intermediate = scope["intermediate"]

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
    tmpdir,
    scope,
    indexes_path,
    snapshot,
):
    temp_indexes_path = Path(tmpdir)/"index"
    temp_indexes_path.mkdir()

    for path in indexes_path.iterdir():
        if "reference" in path.name:
            logger.info(path)
            shutil.copyfile(
                path,
                temp_indexes_path/path.name.replace("reference", "isolates")
            )

    scope["isolate_path"] = temp_indexes_path

    map_isolates = await scope.bind(workflow.map_isolates)
    await map_isolates()

    vta_path = scope["intermediate"].isolate_vta_path

    with vta_path.open("r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data, "isolates")


async def test_map_subtraction(
    tmpdir,
    scope,
    fastq_path,
    host_path,
    snapshot,
):
    temp_subtraction_path = Path(tmpdir) / "subtraction"
    temp_subtraction_path.mkdir()

    mapped_fastq_path = temp_subtraction_path/"mapped.fastq"

    shutil.copyfile(fastq_path, mapped_fastq_path)

    scope["subtraction"] = SimpleNamespace(path=host_path)

    scope["intermediate"] = SimpleNamespace(
        isolate_mapped_fastq_path=mapped_fastq_path,
    )

    map_subtractions = await scope.bind(workflow.map_subtractions)
    await map_subtractions()

    lines = sorted(scope["intermediate"].to_subtraction)

    snapshot.assert_match(lines)


async def test_subtract_mapping(tmpdir, vta_path, to_subtraction_path, scope):
    temp_vta_path = Path(tmpdir)/"to_isolates.vta"
    shutil.copyfile(vta_path, temp_vta_path)

    with to_subtraction_path.open("r") as f:
        scope["intermediate"] = SimpleNamespace(
            to_subtraction=json.load(f),
            isolate_vta_path=temp_vta_path,
        )

    scope["isolate_path"] = Path(tmpdir)

    subtract_mapping = await scope.bind(workflow.subtract_mapping)
    await subtract_mapping()

    assert scope["results"]["subtracted_count"] == 4


async def test_pathoscope(
    tmpdir,
    scope,
    vta_path,
    ref_lengths_path,
    snapshot,
):
    temp_vta_path = Path(tmpdir) / "to_isolates.vta"
    with ref_lengths_path.open("r") as f:
        scope["intermediate"] = SimpleNamespace(
            lengths=json.load(f),
            isolate_vta_path=temp_vta_path,
        )

    shutil.copyfile(vta_path, temp_vta_path)

    scope["isolate_path"] = Path(tmpdir)

    reassignment = await scope.bind(workflow.reassignment)
    await reassignment()

    intermediate = scope["intermediate"]

    with intermediate.reassigned_path.open("r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data, 1)

    with intermediate.report_path.open("r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data, 2)

    snapshot.assert_match(scope["results"], 3)
