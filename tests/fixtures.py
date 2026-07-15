from pathlib import Path

import pytest


@pytest.fixture()
def example_path():
    return Path(__file__).parent.parent / "example"


def get_sam_lines():
    with open(Path(__file__).parent.parent / "example/fifty_lines.sam") as handle:
        return handle.read().split("\n")[0:-1]


@pytest.fixture(params=get_sam_lines(), ids=lambda x: x.split("\t")[0])
def sam_line(request):
    return request.param.split("\t")
