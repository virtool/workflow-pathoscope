from . import _pathoscope
from virtool_workflow import fixture


@fixture
def pathoscope():
    return _pathoscope
