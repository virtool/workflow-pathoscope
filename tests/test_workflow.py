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

        shutil.copyfile(fastq_path.with_suffix(".fq.gz"), temp_read_path/"reads_1.fq.gz")

        _scope["reads"] = Reads(
            sample=_scope["sample"],
            quality={},
            path=temp_read_path
        )

        run_in_executor = await _scope.get_or_instantiate("run_in_executor")
        _scope["run_subprocess"] = run_subprocess = await _scope.get_or_instantiate("_run_subprocess")

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
    scope["intermediate"] = SimpleNamespace(to_otus={})

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


async def test_map_isolates(tmpdir, fastq_path, indexes_path, scope, vta_path, snapshot):
    isolates_path = Path(tmpdir)/"isolates"
    isolates_path.mkdir()

    scope["intermediate"] = SimpleNamespace()
    scope["isolate_path"] = isolates_path

    for filename in os.listdir(indexes_path):
        if "reference" in filename:
            shutil.copyfile(
                indexes_path/filename,
                isolates_path/filename.replace("reference", "isolates")
            )

    map_isolates = await scope.bind(workflow.map_isolates)

    await map_isolates()

    intermediate = scope["intermediate"]

    with open(intermediate.vta_path, "r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data, "isolates")


async def test_map_subtraction(tmpdir, snapshot, host_path, fastq_path, scope):
    scope["subtraction"] = Subtraction(
        id="subtraction1",
        name="subtraction1",
        nickname=None,
        count=None,
        deleted=None,
        gc=None,
        is_host=True,
        path=host_path,
    )

    shutil.copyfile(fastq_path, os.path.join(
        mock_job.params["temp_analysis_path"], "mapped.fastq"))

    await virtool.jobs.pathoscope.map_subtraction(mock_job)

    sorted_lines = sorted(mock_job.intermediate["to_subtraction"])

    snapshot.assert_match(sorted_lines, "subtraction")


async def test_subtract_mapping(mock_job, to_subtraction_path, vta_path):
    with open(to_subtraction_path, "r") as handle:
        mock_job.intermediate["to_subtraction"] = json.load(handle)

    shutil.copyfile(vta_path, os.path.join(
        mock_job.params["temp_analysis_path"], "to_isolates.vta"))

    await virtool.jobs.pathoscope.subtract_mapping(mock_job)

    assert mock_job.results["subtracted_count"] == 4


async def test_pathoscope(snapshot, mock_job, ref_lengths_path, vta_path):
    with open(ref_lengths_path, "r") as handle:
        mock_job.intermediate["ref_lengths"] = json.load(handle)

    shutil.copyfile(
        vta_path,
        os.path.join(mock_job.params["temp_analysis_path"], "to_isolates.vta")
    )

    mock_job.params["sequence_otu_map"] = {
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

    await virtool.jobs.pathoscope.pathoscope(mock_job)

    with open(os.path.join(mock_job.params["temp_analysis_path"], "reassigned.vta"), "r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data)

    with open(os.path.join(mock_job.params["temp_analysis_path"], "report.tsv"), "r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data)

    snapshot.assert_match(mock_job.results)

