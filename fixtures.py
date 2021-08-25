from pathlib import Path
from types import SimpleNamespace
from typing import List

from virtool_workflow.analysis.indexes import Index
from virtool_workflow.data_model import Subtraction

import pathoscope
from virtool_workflow import fixture


@fixture
def pathoscope():
    return pathoscope


@fixture
def index(indexes: List[Index]):
    return indexes[0]


@fixture
def subtraction(subtractions: List[Subtraction]):
    return subtractions[0]


@fixture
def intermediate():
    """A namespace for storing intermediate values."""
    return SimpleNamespace(
        to_otus=set(),
    )


@fixture
def isolate_path(work_path: Path):
    path = work_path / "isolates"
    path.mkdir()

    return path


@fixture
def isolate_fasta_path(isolate_path: Path):
    return isolate_path / "isolate_index.fa"


@fixture
def isolate_fastq_path(isolate_path: Path):
    return isolate_path / "isolate_mapped.fq"


@fixture
def isolate_index_path(isolate_path: Path):
    return isolate_path / "isolates"


@fixture
def isolate_vta_path(isolate_path: Path):
    return isolate_path / "to_isolates.vta"


@fixture
def subtraction_vta_path(work_path: Path):
    return work_path / "subtracted.vta"


@fixture
def p_score_cutoff():
    return 0.01
