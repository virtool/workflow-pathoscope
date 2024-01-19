from pathlib import Path
from types import SimpleNamespace

from pyfixtures import fixture


@fixture
def intermediate():
    """A namespace for storing intermediate values."""
    return SimpleNamespace(
        isolate_high_scores={},
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
def isolate_sam_path(isolate_path: Path):
    return isolate_path / "to_isolates.sam"


@fixture
def p_score_cutoff():
    return 0.01


@fixture
def subtracted_sam_path(work_path: Path) -> Path:
    """The path to the SAM file after subtraction reads have been eliminated."""
    return work_path / "subtracted.sam"
