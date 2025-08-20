"""Test the combined eliminate_subtraction_and_filter_fastq function."""

from pathlib import Path
from workflow_pathoscope.rust import eliminate_subtraction_and_filter_fastq

# Test SAM data for isolates - using correct format
isolate_sam_content = """@HD	VN:1.0	SO:unsorted
@SQ	SN:ref1	LN:1000
@SQ	SN:ref2	LN:2000
read1	0	ref1	100	255	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	AS:i:45
read2	0	ref2	200	255	30M	*	0	0	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	*	AS:i:25
read3	0	ref1	300	255	40M	*	0	0	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	*	AS:i:30
"""

# Test SAM data for subtraction (higher scores for some reads)
subtraction_sam_content = """@HD	VN:1.0	SO:unsorted
@SQ	SN:plant_ref	LN:1000
read1	0	plant_ref	100	255	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	AS:i:60
read3	0	plant_ref	300	255	40M	*	0	0	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	*	AS:i:20
"""

# Test FASTQ data
fastq_content = """@read1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@read3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
@read4
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
"""


def test_combined_function(tmp_path):
    """Test the combined eliminate_subtraction_and_filter_fastq function."""
    # Create test files
    isolate_sam_path = tmp_path / "isolate.sam"
    subtraction_sam_path = tmp_path / "subtraction.sam"
    output_sam_path = tmp_path / "output.sam"
    input_fastq_path = tmp_path / "input.fq"
    output_fastq_path = tmp_path / "output.fq"

    with open(isolate_sam_path, "w") as f:
        f.write(isolate_sam_content)

    with open(subtraction_sam_path, "w") as f:
        f.write(subtraction_sam_content)

    with open(input_fastq_path, "w") as f:
        f.write(fastq_content)

    # Call the combined function
    eliminated_count = eliminate_subtraction_and_filter_fastq(
        str(isolate_sam_path),
        str(subtraction_sam_path),
        str(output_sam_path),
        str(input_fastq_path),
        str(output_fastq_path),
    )

    # Verify the function returned the correct count
    # read1 should be eliminated (subtraction score 60+50=110 > isolate score 45+50=95)
    # read3 should NOT be eliminated (subtraction score 20+40=60 < isolate score 30+40=70)
    # read2 has no subtraction alignment, so kept
    # Only read1 should be eliminated
    assert eliminated_count == 1, (
        f"Expected 1 eliminated read, got {eliminated_count}"
    )

    # Verify output FASTQ file
    with open(output_fastq_path, "r") as f:
        output_content = f.read()

    # read1 should be eliminated
    assert "read1" not in output_content, "read1 should be eliminated from FASTQ"
    # read2, read3, read4 should be kept
    assert "read2" in output_content, "read2 should be kept in FASTQ"
    assert "read3" in output_content, "read3 should be kept in FASTQ"
    assert "read4" in output_content, "read4 should be kept in FASTQ"

    # Verify output BAM file exists
    assert output_sam_path.exists(), "Output BAM file should be created"


def test_combined_function_no_eliminations(tmp_path):
    """Test when no reads are eliminated."""
    # Create test files where subtraction scores are lower
    isolate_sam_path = tmp_path / "isolate.sam"
    subtraction_sam_path = tmp_path / "subtraction.sam"
    output_sam_path = tmp_path / "output.sam"
    input_fastq_path = tmp_path / "input.fq"
    output_fastq_path = tmp_path / "output.fq"

    # Lower subtraction scores
    low_subtraction_content = """@HD	VN:1.0	SO:unsorted
@SQ	SN:plant_ref	LN:1000
read1	0	plant_ref	100	255	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	AS:i:10
"""

    with open(isolate_sam_path, "w") as f:
        f.write(isolate_sam_content)

    with open(subtraction_sam_path, "w") as f:
        f.write(low_subtraction_content)

    with open(input_fastq_path, "w") as f:
        f.write(fastq_content)

    # Call the combined function
    eliminated_count = eliminate_subtraction_and_filter_fastq(
        str(isolate_sam_path),
        str(subtraction_sam_path),
        str(output_sam_path),
        str(input_fastq_path),
        str(output_fastq_path),
    )

    # No reads should be eliminated
    assert eliminated_count == 0, (
        f"Expected 0 eliminated reads, got {eliminated_count}"
    )

    # All reads should be in output FASTQ
    with open(output_fastq_path, "r") as f:
        output_content = f.read()

    assert "read1" in output_content, "read1 should be kept"
    assert "read2" in output_content, "read2 should be kept"
    assert "read3" in output_content, "read3 should be kept"
    assert "read4" in output_content, "read4 should be kept"
