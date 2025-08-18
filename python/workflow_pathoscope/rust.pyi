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
    sam_path: str,
    p_score_cutoff: float,
    ref_lengths: dict[str, int],
) -> PathoscopeResults:
    """Run Pathoscope expectation maximization algorithm using Rust."""
