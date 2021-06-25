import pytest


@pytest.fixture
def best_hit_path(test_files_path):
    return test_files_path / "best_hit"


@pytest.fixture
def results_path(test_files_path):
    return test_files_path / "results.json"


@pytest.fixture
def em_path(test_files_path):
    return test_files_path / "em"


@pytest.fixture
def isolates_vta_path(test_files_path):
    return test_files_path / "to_isolates.vta"


@pytest.fixture
def matrix_path(test_files_path):
    return test_files_path / "ps_matrix"


@pytest.fixture
def ref_lengths_path(test_files_path):
    return test_files_path / "ref_lengths.json"


@pytest.fixture
def sam_path(test_files_path):
    return test_files_path / "test_al.sam"


@pytest.fixture
def scores(test_files_path):
    return test_files_path / "scores"


@pytest.fixture
def to_subtraction_path(test_files_path):
    return test_files_path / "to_subtraction.json"


@pytest.fixture
def unu_path(test_files_path):
    return test_files_path / "unu"


@pytest.fixture
def vta_path(test_files_path):
    return test_files_path / "test.vta"


@pytest.fixture
def indexes_path(test_files_path):
    return test_files_path / "index"


@pytest.fixture
def fastq_path(test_files_path):
    return test_files_path / "test.fq"


@pytest.fixture
def host_path(indexes_path):
    return indexes_path / "host"
