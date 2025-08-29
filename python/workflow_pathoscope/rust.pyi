def run_eliminate_subtraction(
    isolate_sam_path: str, subtraction_sam_path: str, output_sam_path: str
) -> None:
    """Eliminate subtraction reads from isolate reads using Rust."""

class PathoscopeResults:
    best_hit_initial_reads: list[float]
    best_hit_initial: list[float]
    level_1_initial: list[float]
    level_2_initial: list[float]
    best_hit_final_reads: list[float]
    best_hit_final: list[float]
    level_1_final: list[float]
    level_2_final: list[float]
    init_pi: list[float]
    pi: list[float]
    refs: list[str]
    reads: list[str]
    coverage: dict[str, list[int]]

def run_expectation_maximization(
    alignment_path: str,
    p_score_cutoff: float,
    ref_lengths: dict[str, int],
) -> PathoscopeResults:
    """Run Pathoscope expectation maximization algorithm using Rust on SAM/BAM files."""

def parse_isolate_scores(
    alignment_path: str,
    p_score_cutoff: float,
) -> dict[str, float]:
    """Parse isolate alignment file (SAM or BAM) and extract high scores for each read."""

def find_candidate_otus(
    alignment_path: str,
    p_score_cutoff: float,
) -> set[str]:
    """Extract candidate OTU reference IDs from an alignment file (SAM/BAM)."""

def find_candidate_otus_from_bytes(
    sam_bytes: bytes,
    p_score_cutoff: float,
) -> set[str]:
    """Extract candidate OTU reference IDs from SAM text data."""

def calculate_coverage_from_em_results(
    alignment_path: str,
    p_score_cutoff: float,
    ref_lengths: dict[str, int],
) -> dict[str, list[int]]:
    """Calculate coverage directly from EM results and alignment data."""
