import shutil
import pytest
from pathlib import Path


@pytest.fixture
def test_files_path():
    return Path(__file__).parent / "test_files"


@pytest.fixture
def test_sam_path(test_files_path, tmp_path):
    src_path = test_files_path / "test_al.sam"
    dst_path = tmp_path / "test_sam_file"
    dst_path.mkdir()
    dst_path = dst_path / "test_al.sam"
    shutil.copy(src_path, dst_path)

    return dst_path


def get_sam_lines():
    path = Path(__file__).parent / "test_files" / "sam_50.sam"

    with open(path, "r") as handle:
        return handle.read().split("\n")[0:-1]


@pytest.fixture(params=get_sam_lines(), ids=lambda x: x.split("\t")[0])
def sam_line(request):
    return request.param.split("\t")
