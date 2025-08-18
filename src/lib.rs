mod subtraction;
mod coverage;
mod matrix;
mod em;
mod parse_sam;

use subtraction::eliminate_subtraction;
use matrix::build_matrix;
use em::{em, compute_best_hit};
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use std::collections::HashMap;

// Type aliases for complex HashMap types used throughout the codebase
pub type UniqueReads = HashMap<i32, (i32, f64)>;
pub type MultiMappingReads = HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;

#[pyclass]
#[derive(Clone)]
pub struct PathoscopeResults {
    #[pyo3(get)]
    pub best_hit_initial_reads: Vec<f64>,
    #[pyo3(get)]
    pub best_hit_initial: Vec<f64>,
    #[pyo3(get)]
    pub level_1_initial: Vec<f64>,
    #[pyo3(get)]
    pub level_2_initial: Vec<f64>,
    #[pyo3(get)]
    pub best_hit_final_reads: Vec<f64>,
    #[pyo3(get)]
    pub best_hit_final: Vec<f64>,
    #[pyo3(get)]
    pub level_1_final: Vec<f64>,
    #[pyo3(get)]
    pub level_2_final: Vec<f64>,
    #[pyo3(get)]
    pub init_pi: Vec<f64>,
    #[pyo3(get)]
    pub pi: Vec<f64>,
    #[pyo3(get)]
    pub refs: Vec<String>,
    #[pyo3(get)]
    pub reads: Vec<String>,
    #[pyo3(get)]
    pub coverage: HashMap<String, Vec<usize>>,
}

#[pymodule]
/// pyo3 interface
fn rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PathoscopeResults>()?;
    m.add_function(wrap_pyfunction!(run_expectation_maximization, m)?)?;
    m.add_function(wrap_pyfunction!(run_eliminate_subtraction, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_coverage_from_em_results, m)?)?;
    Ok(())
}

#[pyfunction]
/// Calculate coverage directly from EM results and SAM data
pub fn calculate_coverage_from_em_results(
    _py: Python,
    sam_path: String,
    p_score_cutoff: f64,
    ref_lengths: HashMap<String, usize>,
) -> PyResult<HashMap<String, Vec<usize>>> {
    // Build matrix to get EM input data and minimal alignments
    let (u, nu, _refs, _reads, minimal_alignments) =
        build_matrix(sam_path.as_str(), None)
            .map_err(|e| PyErr::new::<PyIOError, _>(format!("Failed to build matrix: {}", e)))?;

    // Run EM algorithm
    let (_init_pi, _pi, _theta, nu) = em(&u, nu, &_refs, 50, 1e-7, 0.0, 0.0);

    // Calculate coverage directly
    let coverage = coverage::calculate_coverage_from_em(
        &u,
        &nu,
        &minimal_alignments,
        &_refs,
        p_score_cutoff,
        &ref_lengths,
    );

    Ok(coverage)
}

#[pyfunction]
/// Entry point for eliminate_subtraction - PyO3 wrapper
pub fn run_eliminate_subtraction(
    _py: Python,
    isolate_sam_path: String,
    subtraction_sam_path: String,
    output_sam_path: String,
) -> PyResult<()> {
    // Call the pure Rust function and map errors to PyResult
    eliminate_subtraction(&isolate_sam_path, &subtraction_sam_path, &output_sam_path)
        .map_err(|e| PyErr::new::<PyIOError, _>(e.to_string()))
}


#[pyfunction]
/// Entry point for expectation_maximization
pub fn run_expectation_maximization(
    _py: Python,
    sam_path: String,
    p_score_cutoff: f64,
    ref_lengths: HashMap<String, usize>,
) -> PyResult<PathoscopeResults> {
    let (u, nu, refs, reads, minimal_alignments) =
        build_matrix(sam_path.as_str(), None)
            .map_err(|e| PyErr::new::<PyIOError, _>(format!("Failed to build matrix: {}", e)))?;

    let (best_hit_initial_reads, best_hit_initial, level_1_initial, level_2_initial) =
        compute_best_hit(&u, &nu, &refs, &reads);

    let (init_pi, pi, _, nu) = em(&u, nu, &refs, 50, 1e-7, 0.0, 0.0);

    let (best_hit_final_reads, best_hit_final, level_1_final, level_2_final) =
        compute_best_hit(&u, &nu, &refs, &reads);

    // Calculate coverage directly from EM results
    let coverage = coverage::calculate_coverage_from_em(
        &u,
        &nu,
        &minimal_alignments,
        &refs,
        p_score_cutoff,
        &ref_lengths,
    );

    Ok(PathoscopeResults {
        best_hit_initial_reads,
        best_hit_initial,
        level_1_initial,
        level_2_initial,
        best_hit_final_reads,
        best_hit_final,
        level_1_final,
        level_2_final,
        init_pi,
        pi,
        refs,
        reads,
        coverage,
    })
}








/// Tests
#[cfg(test)]
mod tests {
    #![allow(unused)]

    use crate::*;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    use std::io::Read;

    extern crate yaml_rust;
    use yaml_rust::{YamlEmitter, YamlLoader};

    #[test]
    fn test_compute_best_hit() {
        let mut unique_reads: UniqueReads = std::collections::HashMap::new();
        let mut multi_reads: MultiMappingReads = std::collections::HashMap::new();
        
        // Unique reads: each maps to exactly one reference
        let read0_maps_to_ref0 = (0, 100.0);
        let read1_maps_to_ref1 = (1, 100.0);
        unique_reads.insert(0, read0_maps_to_ref0);
        unique_reads.insert(1, read1_maps_to_ref1);
        
        // Multi-mapping reads: each maps to multiple references with different scores
        let read2_prefers_ref0 = (vec![0, 1], vec![90.0, 10.0], vec![0.9, 0.1], 90.0);
        let read3_tied_between_refs = (vec![0, 1], vec![80.0, 80.0], vec![0.5, 0.5], 80.0);
        let read4_tests_thresholds = (vec![0, 1, 2], vec![60.0, 35.0, 5.0], vec![0.6, 0.35, 0.05], 60.0);
        let read5_below_threshold = (vec![0, 1], vec![99.5, 0.5], vec![0.995, 0.005], 99.5);
        
        multi_reads.insert(2, read2_prefers_ref0);
        multi_reads.insert(3, read3_tied_between_refs);
        multi_reads.insert(4, read4_tests_thresholds);
        multi_reads.insert(5, read5_below_threshold);
        
        let refs = vec!["ref0".to_string(), "ref1".to_string(), "ref2".to_string()];
        let reads = vec![
            "read0".to_string(), "read1".to_string(), "read2".to_string(), 
            "read3".to_string(), "read4".to_string(), "read5".to_string()
        ];
        
        let (best_hit_reads, best_hit, level1, level2) = em::compute_best_hit(&unique_reads, &multi_reads, &refs, &reads);
        
        // Test 1: Verify best hit assignments (raw counts)
        let expected_ref0_best_hits = 1.0   // read0 (unique)
                                    + 1.0   // read2 (0.9 > 0.1)
                                    + 0.5   // read3 (tied 0.5 == 0.5)
                                    + 1.0   // read4 (0.6 > 0.35 > 0.05)
                                    + 1.0;  // read5 (0.995 > 0.005)
        assert!((best_hit_reads[0] - expected_ref0_best_hits).abs() < 0.001, 
                "ref0 best_hit_reads should be {}, got {}", expected_ref0_best_hits, best_hit_reads[0]);
        
        let expected_ref1_best_hits = 1.0   // read1 (unique)
                                    + 0.0   // read2 (0.1 < 0.9)
                                    + 0.5   // read3 (tied 0.5 == 0.5)
                                    + 0.0   // read4 (0.35 < 0.6)
                                    + 0.0;  // read5 (0.005 < 0.995)
        assert!((best_hit_reads[1] - expected_ref1_best_hits).abs() < 0.001,
                "ref1 best_hit_reads should be {}, got {}", expected_ref1_best_hits, best_hit_reads[1]);
        
        let expected_ref2_best_hits = 0.0;  // No reads have ref2 as best hit
        assert!((best_hit_reads[2] - expected_ref2_best_hits).abs() < 0.001,
                "ref2 best_hit_reads should be {}, got {}", expected_ref2_best_hits, best_hit_reads[2]);
        
        // Test 2: Verify normalization (divided by total read count = 6)
        let total_reads = 6.0;
        assert!((best_hit[0] - expected_ref0_best_hits/total_reads).abs() < 0.001,
                "ref0 normalized best_hit incorrect");
        assert!((best_hit[1] - expected_ref1_best_hits/total_reads).abs() < 0.001,
                "ref1 normalized best_hit incorrect");
        assert!((best_hit[2] - expected_ref2_best_hits/total_reads).abs() < 0.001,
                "ref2 normalized best_hit incorrect");
        
        // Test 3: Verify level1 assignments (score >= 0.5)
        let expected_ref0_level1 = 1.0   // read0 (unique, always counts)
                                 + 1.0   // read2 (0.9 >= 0.5)
                                 + 1.0   // read3 (0.5 >= 0.5)
                                 + 1.0   // read4 (0.6 >= 0.5)
                                 + 1.0;  // read5 (0.995 >= 0.5)
        assert!((level1[0] - expected_ref0_level1/total_reads).abs() < 0.001,
                "ref0 level1 should be {}", expected_ref0_level1/total_reads);
        
        let expected_ref1_level1 = 1.0   // read1 (unique, always counts)
                                 + 0.0   // read2 (0.1 < 0.5)
                                 + 1.0   // read3 (0.5 >= 0.5)
                                 + 0.0   // read4 (0.35 < 0.5)
                                 + 0.0;  // read5 (0.005 < 0.5)
        assert!((level1[1] - expected_ref1_level1/total_reads).abs() < 0.001,
                "ref1 level1 should be {}", expected_ref1_level1/total_reads);
        
        let expected_ref2_level1 = 0.0;  // read4 (0.05 < 0.5)
        assert!((level1[2] - expected_ref2_level1/total_reads).abs() < 0.001,
                "ref2 level1 should be 0.0");
        
        // Test 4: Verify level2 assignments (0.01 <= score < 0.5)
        let expected_ref0_level2 = 0.0;  // All ref0 scores >= 0.5, so counted in level1
        assert!((level2[0] - expected_ref0_level2/total_reads).abs() < 0.001,
                "ref0 level2 should be 0.0");
        
        let expected_ref1_level2 = 1.0   // read2 (0.1 >= 0.01 and < 0.5)
                                 + 0.0   // read3 (0.5 >= 0.5, counted in level1)
                                 + 1.0   // read4 (0.35 >= 0.01 and < 0.5)
                                 + 0.0;  // read5 (0.005 < 0.01)
        assert!((level2[1] - expected_ref1_level2/total_reads).abs() < 0.001,
                "ref1 level2 should be {}", expected_ref1_level2/total_reads);
        
        let expected_ref2_level2 = 1.0;  // read4 (0.05 >= 0.01 and < 0.5)
        assert!((level2[2] - expected_ref2_level2/total_reads).abs() < 0.001,
                "ref2 level2 should be {}", expected_ref2_level2/total_reads);
        
        // Test 5: Verify output vector lengths
        assert_eq!(best_hit_reads.len(), 3, "best_hit_reads should have 3 elements");
        assert_eq!(best_hit.len(), 3, "best_hit should have 3 elements");
        assert_eq!(level1.len(), 3, "level1 should have 3 elements");
        assert_eq!(level2.len(), 3, "level2 should have 3 elements");
    }

    #[test]
    fn test_compute_best_hit_edge_cases() {
        // Test empty inputs
        let u: UniqueReads = std::collections::HashMap::new();
        let nu: MultiMappingReads = std::collections::HashMap::new();
        let refs = vec!["ref0".to_string()];
        let reads = vec!["read0".to_string()];
        
        let (best_hit_reads, best_hit, level1, level2) = em::compute_best_hit(&u, &nu, &refs, &reads);
        
        assert_eq!(best_hit_reads[0], 0.0, "Empty input should result in zero counts");
        assert_eq!(best_hit[0], 0.0, "Empty input should result in zero normalized values");
        assert_eq!(level1[0], 0.0, "Empty input should result in zero level1");
        assert_eq!(level2[0], 0.0, "Empty input should result in zero level2");
        
        // Test with nu entry that has no valid scores
        let mut nu_invalid: MultiMappingReads = std::collections::HashMap::new();
        nu_invalid.insert(0, (vec![0], vec![10.0], vec![0.0], 10.0)); // normalized score is 0.0
        
        let reads_single = vec!["read0".to_string()];
        let (best_hit_reads, best_hit, level1, level2) = em::compute_best_hit(&u, &nu_invalid, &refs, &reads_single);
        
        // Even with 0.0 normalized score, it should still be considered as best hit (0.0 == max)
        assert_eq!(best_hit_reads[0], 1.0, "Zero score should still count as best hit when it's the only/max score");
    }

    #[test]
    fn test_eliminate_subtraction_integration() {
        use std::fs;
        use tempfile::TempDir;

        // Create temp directory for output
        let temp_dir = TempDir::new().unwrap();
        let output_sam_path = temp_dir.path().join("output.sam");

        // Run the pure Rust function directly
        subtraction::eliminate_subtraction(
            "example/rust/test_isolates_minimal.sam",
            "example/rust/test_subtraction_minimal.sam",
            output_sam_path.to_str().unwrap(),
        )
        .unwrap();

        // Read and verify output SAM
        let output_content = fs::read_to_string(&output_sam_path).unwrap();
        let output_lines: Vec<&str> = output_content.lines().collect();

        // Verify headers are preserved
        assert!(output_lines[0].starts_with("@HD"), "Headers should be preserved");
        assert!(
            output_lines[1].starts_with("@SQ\tSN:ref1"),
            "First sequence header should be preserved"
        );
        assert!(
            output_lines[2].starts_with("@SQ\tSN:ref2"),
            "Second sequence header should be preserved"
        );
        assert!(
            output_lines[3].starts_with("@PG"),
            "Program header should be preserved"
        );
        assert_eq!(
            output_lines[4], "@CO\tTest comment line",
            "Comment lines should be preserved in proper SAM format"
        );

        // Verify correct reads are kept
        assert!(
            output_content.contains("read_keep"),
            "Read with lower subtraction score should be kept"
        );
        assert!(
            output_content.contains("read_unknown"),
            "Read not in subtraction should be kept"
        );

        // Verify eliminated reads are not present
        assert!(
            !output_content.contains("read_eliminate"),
            "Read with higher subtraction score should be eliminated"
        );
        assert!(
            !output_content.contains("read_equal"),
            "Read with equal scores should be eliminated"
        );
        assert!(
            !output_content.contains("read_unmapped"),
            "Unmapped reads should be skipped"
        );

        // Verify subtracted_read_ids.txt
        let ids_path = temp_dir.path().join("subtracted_read_ids.txt");
        let ids_content = fs::read_to_string(&ids_path).unwrap();
        assert!(
            ids_content.contains("read_eliminate"),
            "Eliminated read should be in subtracted_read_ids.txt"
        );
        assert!(
            ids_content.contains("read_equal"),
            "Equal score read should be in subtracted_read_ids.txt"
        );
        assert!(
            !ids_content.contains("read_keep"),
            "Kept read should not be in subtracted_read_ids.txt"
        );
        assert!(
            !ids_content.contains("read_unknown"),
            "Unknown read should not be in subtracted_read_ids.txt"
        );
    }
}
