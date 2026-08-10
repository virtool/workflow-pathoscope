import json
from pathlib import Path
from types import SimpleNamespace

from syrupy import SnapshotAssertion
from workflow_pathoscope.utils import (
    write_report,
)


def test_write_report(
    snapshot: SnapshotAssertion,
    tmp_path,
):
    """Test that a report is written correctly given a set of Pathoscope results."""
    report_path = tmp_path / "report.tsv"

    with open(Path(__file__).parent / "report_input.json") as f:
        data = json.load(f)

    # Create a mock PathoscopeResults object from the test data
    pathoscope_results = SimpleNamespace(
        best_hit_initial=data["best_hit_initial"],
        best_hit_initial_reads=data["best_hit_initial_reads"],
        best_hit_final=data["best_hit_final"],
        best_hit_final_reads=data["best_hit_final_reads"],
        init_pi=data["init_pi"],
        pi=data["pi"],
        refs=data["refs"],
        reads=["dummy_read"] * data["read_count"],  # Create dummy reads list
        level_1_initial=data["level_1_initial"],
        level_2_initial=data["level_2_initial"],
        level_1_final=data["level_1_final"],
        level_2_final=data["level_2_final"],
        coverage={},  # Not used in write_report
    )

    write_report(report_path, pathoscope_results)

    with open(report_path) as f:
        assert f.read() == snapshot
