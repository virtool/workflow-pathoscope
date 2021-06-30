import logging
import json
import os
import shutil

import pytest

from virtool_workflow_pathoscope import workflow, fixtures
from virtool_workflow import FixtureScope
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.analysis.library_types import LibraryType
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.subtractions import Subtraction
from virtool_workflow.data_model.samples import Sample
from pathlib import Path
from types import SimpleNamespace
from virtool_workflow.runtime.fixtures import analysis as analysis_fixtures


logger = logging.getLogger(__name__)


@pytest.fixture
async def scope(tmpdir, fastq_path, otu_resource, indexes_path):
    """Fixture scope for pathoscope workflow tests"""
    async with FixtureScope(analysis_fixtures) as _scope:
        _scope.fixture(fixtures.pathoscope)
        _scope.fixture(workflow.p_score_cutoff)
        _scope.fixture(workflow.intermediate)

        # Create a mock sample fixture
        _scope["sample"] = Sample(
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

        _scope["reads"] = Reads(
            sample=_scope["sample"],
            quality={},
            path=temp_read_path
        )

        run_in_executor = await _scope.get_or_instantiate("run_in_executor")
        run_subprocess = await _scope.get_or_instantiate("_run_subprocess")
        _scope["run_subprocess"] = run_subprocess

        # Create a mock index fixture
        temp_index_path = Path(tmpdir) / "indexes"
        shutil.copytree(indexes_path, temp_index_path)

        _scope["index"] = Index(
            id="index3",
            manifest={
                "foobar": 10,
                "reo": 5,
                "baz": 6,
            },
            _sequence_otu_map=otu_resource[0],
            reference=None,
            path=temp_index_path,
            ready=True,
            _run_in_executor=run_in_executor,
            _run_subprocess=run_subprocess,
        )

        _scope["logger"] = logger

        yield _scope


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


async def test_map_isolates(tmpdir, scope, indexes_path, snapshot):
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


async def test_map_subtraction(tmpdir, scope, fastq_path, host_path, snapshot):
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
