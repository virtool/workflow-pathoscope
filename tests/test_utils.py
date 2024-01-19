import json
import shutil
from pathlib import Path

from syrupy import SnapshotAssertion

from utils import calculate_coverage, find_sam_align_score, write_report


def test_write_report(
    snapshot: SnapshotAssertion,
    tmp_path,
):
    """Test that a report is written correctly given a set of Pathoscope results."""
    report_path = tmp_path / "report.tsv"

    with open(Path(__file__).parent / "test_utils" / "report_input.json") as f:
        data = json.load(f)

    write_report(report_path, **data)

    assert open(report_path).read() == snapshot


def test_calculate_coverage(example_path: Path, tmp_path):
    ref_lengths = dict()

    sam_path = tmp_path / "mapped.sam"

    shutil.copyfile(example_path / "to_isolates.sam", sam_path)

    with open(sam_path) as handle:
        for line in handle:
            if line[0:3] == "@SQ":
                ref_id = None
                length = None

                for field in line.rstrip().split("\t"):
                    if "SN:" in field:
                        ref_id = field.split(":")[1]
                    if "LN:" in field:
                        length = int(field.split(":")[1])

                assert ref_id and length
                assert ref_id not in ref_lengths

                ref_lengths[ref_id] = length

    calculate_coverage(sam_path, ref_lengths)


def test_find_sam_align_score(sam_line, snapshot: SnapshotAssertion):
    assert find_sam_align_score(sam_line) == snapshot
