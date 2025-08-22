"""Test the new Rust subtract_fastq function."""

import tempfile
from pathlib import Path
from workflow_pathoscope.rust import subtract_fastq

# Test FASTQ content
test_fastq_content = """@read1
ACGTACGTACGT
+
IIIIIIIIIIII
@read2
TGCATGCATGCA 
+
JJJJJJJJJJJJ
@read3
AAAAAAAAAAAAA
+
KKKKKKKKKKKKK
@read4_with_spaces description here
CCCCCCCCCCCCC
+
LLLLLLLLLLLLL
"""


def test_rust_subtract_fastq():
    """Test the Rust subtract_fastq function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create input FASTQ file
        input_path = Path(tmpdir) / "input.fq"
        output_path = Path(tmpdir) / "output.fq"

        with open(input_path, "w") as f:
            f.write(test_fastq_content)

        # Subtract read2 and read4_with_spaces
        subtracted_reads = {"read2", "read4_with_spaces"}

        # Call Rust function
        kept_count = subtract_fastq(str(input_path), str(output_path), subtracted_reads)

        # Read output file
        with open(output_path, "r") as f:
            output_content = f.read()

        # Verify results
        assert kept_count == 2, f"Expected 2 kept reads, got {kept_count}"
        assert "read1" in output_content, "read1 should be kept"
        assert "read3" in output_content, "read3 should be kept"
        assert "read2" not in output_content, "read2 should be subtracted"
        assert "read4_with_spaces" not in output_content, (
            "read4_with_spaces should be subtracted"
        )


def test_rust_subtract_fastq_empty_subtraction():
    """Test with no reads to subtract."""
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "input.fq"
        output_path = Path(tmpdir) / "output.fq"

        with open(input_path, "w") as f:
            f.write(test_fastq_content)

        # No reads to subtract
        subtracted_reads = set()

        kept_count = subtract_fastq(str(input_path), str(output_path), subtracted_reads)

        with open(output_path, "r") as f:
            output_content = f.read()

        # All reads should be kept
        assert kept_count == 4, f"Expected 4 kept reads, got {kept_count}"
        assert "read1" in output_content
        assert "read2" in output_content
        assert "read3" in output_content
        assert "read4_with_spaces" in output_content


def test_rust_subtract_fastq_all_subtracted():
    """Test subtracting all reads."""
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "input.fq"
        output_path = Path(tmpdir) / "output.fq"

        with open(input_path, "w") as f:
            f.write(test_fastq_content)

        # Subtract all reads
        subtracted_reads = {"read1", "read2", "read3", "read4_with_spaces"}

        kept_count = subtract_fastq(str(input_path), str(output_path), subtracted_reads)

        with open(output_path, "r") as f:
            output_content = f.read()

        # No reads should be kept
        assert kept_count == 0, f"Expected 0 kept reads, got {kept_count}"
        assert output_content == "", "Output file should be empty"
