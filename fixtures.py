from pathlib import Path
from types import SimpleNamespace

from pyfixtures import fixture
from virtool.workflow.data.indexes import INDEX_SQLITE_FILE_NAME


def get_reference_index_path(work_path: Path) -> Path:
    return work_path / "reference_index" / "reference"


def get_collapsed_reference_path(work_path: Path) -> Path:
    return work_path / "collapsed_reference" / INDEX_SQLITE_FILE_NAME


def get_subtraction_indexes_path(work_path: Path) -> Path:
    return work_path / "subtraction_indexes"


def get_isolate_path(work_path: Path) -> Path:
    return work_path / "isolates"


def get_isolate_fastq_path(isolate_path: Path) -> Path:
    return isolate_path / "isolate_mapped.fq"


def get_isolate_index_path(isolate_path: Path) -> Path:
    return isolate_path / "isolates"


def get_isolate_bam_path(isolate_path: Path) -> Path:
    return isolate_path / "to_isolates.bam"


def get_subtracted_bam_path(work_path: Path) -> Path:
    return work_path / "subtracted.bam"


@fixture
def intermediate():
    """A namespace for storing intermediate values."""
    return SimpleNamespace(
        candidate_sequence_ids=set(),
    )


@fixture
def isolate_path(work_path: Path):
    path = get_isolate_path(work_path)
    path.mkdir()

    return path


@fixture
def reference_index_path(work_path: Path):
    return get_reference_index_path(work_path)


@fixture
def collapsed_reference_path(work_path: Path):
    return get_collapsed_reference_path(work_path)


@fixture
def subtraction_indexes_path(work_path: Path) -> Path:
    return get_subtraction_indexes_path(work_path)


@fixture
def isolate_fasta_path(isolate_path: Path):
    return isolate_path / "isolate_index.fa"


@fixture
def isolate_fastq_path(isolate_path: Path):
    return get_isolate_fastq_path(isolate_path)


@fixture
def isolate_index_path(isolate_path: Path):
    return get_isolate_index_path(isolate_path)


@fixture
def isolate_bam_path(isolate_path: Path):
    return get_isolate_bam_path(isolate_path)


@fixture
def p_score_cutoff():
    return 0.01


@fixture
def subtracted_bam_path(work_path: Path) -> Path:
    """The path to the BAM file after subtraction reads have been eliminated."""
    return get_subtracted_bam_path(work_path)
