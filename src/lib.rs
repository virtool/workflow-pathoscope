mod candidates;
mod coverage;
mod em;
mod logging;
mod matrix;
mod sam;
mod subtraction;

use em::{compute_best_hit, em};
use log::info;
use logging::init_logging;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashSet;

use crate::coverage::calculate_coverage_from_bam;

// Type aliases for complex HashMap types used throughout the codebase
pub type UniqueReads = FxHashMap<i32, (i32, f64)>;
pub type MultiMappingReads = FxHashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;
pub type MatrixResult = (
    UniqueReads,
    MultiMappingReads,
    Vec<String>,
    Vec<String>,
    Vec<sam::MinimalAlignment>,
);

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
    pub coverage: FxHashMap<String, Vec<usize>>,
}

/// Lightweight EM results structure that avoids storing full read/ref vectors
/// until absolutely necessary for the final Python interface
#[derive(Clone)]
pub struct CompactEMResults {
    pub best_hit_initial_reads: Vec<f64>,
    pub best_hit_initial: Vec<f64>,
    pub level_1_initial: Vec<f64>,
    pub level_2_initial: Vec<f64>,
    pub best_hit_final_reads: Vec<f64>,
    pub best_hit_final: Vec<f64>,
    pub level_1_final: Vec<f64>,
    pub level_2_final: Vec<f64>,
    pub init_pi: Vec<f64>,
    pub pi: Vec<f64>,
    // Store only reference count and read count to save memory
    pub ref_count: usize,
    pub read_count: usize,
}

impl CompactEMResults {
    /// Convert to full PathoscopeResults when needed
    pub fn to_pathoscope_results(
        self,
        refs: Vec<String>,
        reads: Vec<String>,
        coverage: FxHashMap<String, Vec<usize>>,
    ) -> PathoscopeResults {
        PathoscopeResults {
            best_hit_initial_reads: self.best_hit_initial_reads,
            best_hit_initial: self.best_hit_initial,
            level_1_initial: self.level_1_initial,
            level_2_initial: self.level_2_initial,
            best_hit_final_reads: self.best_hit_final_reads,
            best_hit_final: self.best_hit_final,
            level_1_final: self.level_1_final,
            level_2_final: self.level_2_final,
            init_pi: self.init_pi,
            pi: self.pi,
            refs,
            reads,
            coverage,
        }
    }
}

#[pymodule]
/// pyo3 interface
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PathoscopeResults>()?;
    m.add_function(wrap_pyfunction!(init_logging, m)?)?;
    m.add_function(wrap_pyfunction!(run_expectation_maximization, m)?)?;
    m.add_function(wrap_pyfunction!(find_candidate_otus_with_bowtie2, m)?)?;
    m.add_function(wrap_pyfunction!(run_eliminate_subtraction, m)?)?;
    Ok(())
}

#[pyfunction]
/// Run expectation maximization algorithm
///
/// # Arguments
/// * `alignment_path` - Path to the SAM/BAM file
/// * `p_score_cutoff` - Minimum score threshold for alignments
pub fn run_expectation_maximization(
    py: Python,
    alignment_path: String,
    p_score_cutoff: f64,
) -> PyResult<PathoscopeResults> {
    info!(
        "starting em algorithm: file={}, cutoff={}",
        alignment_path, p_score_cutoff
    );

    py.allow_threads(|| {
        let matrix =
            matrix::build_matrix(alignment_path.as_str(), Some(p_score_cutoff))
                .map_err(|e| {
                    PyErr::new::<PyIOError, _>(format!("Failed to build matrix: {}", e))
                })?;

        // Calculate initial best hit statistics using the matrix
        let initial_best_hit = compute_best_hit(&matrix);

        // Run EM algorithm using the matrix
        let em_results = em(&matrix, 50, 1e-7, 0.0, 0.0);

        // Calculate final best hit statistics using the updated matrix
        let final_best_hit = compute_best_hit(&em_results.updated_matrix);

        // Calculate coverage using the matrix
        let coverage = calculate_coverage_from_bam(
            &alignment_path,
            &em_results.updated_matrix,
            p_score_cutoff,
        )
        .map_err(|e| {
            PyErr::new::<PyIOError, _>(format!("Failed to calculate coverage: {}", e))
        })?;

        Ok(PathoscopeResults {
            best_hit_initial_reads: initial_best_hit.best_hit_reads,
            best_hit_initial: initial_best_hit.best_hit,
            level_1_initial: initial_best_hit.level1,
            level_2_initial: initial_best_hit.level2,
            best_hit_final_reads: final_best_hit.best_hit_reads,
            best_hit_final: final_best_hit.best_hit,
            level_1_final: final_best_hit.level1,
            level_2_final: final_best_hit.level2,
            init_pi: em_results.init_pi,
            pi: em_results.pi,
            refs: matrix.refs,
            reads: matrix.reads,
            coverage,
        })
    })
}

#[pyfunction]
/// Extract candidate OTU reference IDs by running bowtie2 directly with streaming
///
/// This is a PyO3 wrapper around the function in the candidates module.
pub fn find_candidate_otus_with_bowtie2(
    py: Python,
    bowtie_index_path: &str,
    read_paths: Vec<String>,
    proc: i32,
    p_score_cutoff: f64,
) -> PyResult<HashSet<String>> {
    candidates::find_candidate_otus_with_bowtie2(
        py,
        bowtie_index_path,
        read_paths,
        proc,
        p_score_cutoff,
    )
}

#[pyfunction]
/// Eliminate subtraction reads from BAM and filter the input FASTQ file.
///
/// # Arguments
/// * `isolate_sam_path` - Path to the isolate SAM/BAM file
/// * `subtraction_sam_path` - Path to the subtraction SAM file
/// * `output_sam_path` - Path to write the filtered SAM/BAM file
/// * `input_fastq_path` - Path to the input FASTQ file to filter
/// * `output_fastq_path` - Path to write the filtered FASTQ file
/// * `proc` - Number of threads to use for processing
///
/// # Returns
/// Number of reads that were subtracted (eliminated)
pub fn run_eliminate_subtraction(
    py: Python,
    isolate_sam_path: String,
    subtraction_sam_path: String,
    output_sam_path: String,
    input_fastq_path: String,
    output_fastq_path: String,
    proc: i32,
) -> PyResult<usize> {
    // Release the GIL during the CPU-intensive file I/O and processing
    py.allow_threads(|| {
        let proc_usize = proc.max(1) as usize;
        subtraction::eliminate_subtraction(
            &isolate_sam_path,
            &subtraction_sam_path,
            &output_sam_path,
            &input_fastq_path,
            &output_fastq_path,
            proc_usize,
        )
        .map_err(|e| PyErr::new::<PyIOError, _>(e.to_string()))
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
    fn test_run_eliminate_subtraction() {
        use rust_htslib::{bam, bam::Read};
        use std::fs;
        use tempfile::TempDir;

        // Create temp directory for output
        let temp_dir = TempDir::new().unwrap();
        let output_bam_path = temp_dir.path().join("output.bam");
        let input_fastq_path = temp_dir.path().join("input.fastq");
        let output_fastq_path = temp_dir.path().join("output.fastq");

        // Create a minimal FASTQ file for testing
        std::fs::write(&input_fastq_path, "@read_keep\nACGT\n+\nIIII\n@read_eliminate\nACGT\n+\nIIII\n@read_unknown\nACGT\n+\nIIII\n@read_equal\nACGT\n+\nIIII\n").unwrap();

        // Run the pure Rust function directly
        let subtracted_count = subtraction::eliminate_subtraction(
            "example/rust/test_isolates_minimal.sam",
            "example/rust/test_subtraction_minimal.sam",
            output_bam_path.to_str().unwrap(),
            input_fastq_path.to_str().unwrap(),
            output_fastq_path.to_str().unwrap(),
            4,
        )
        .unwrap();

        // Read and verify output BAM
        let mut bam_reader = bam::Reader::from_path(&output_bam_path).unwrap();
        let header = bam_reader.header();

        // Verify headers are preserved
        assert!(
            header.target_count() > 0,
            "Sequence headers should be preserved"
        );

        // Collect read names from BAM output
        let mut read_names = Vec::with_capacity(10);
        let mut record = bam::Record::new();
        while let Some(result) = bam_reader.read(&mut record) {
            if result.is_ok() {
                if let Ok(qname) = std::str::from_utf8(record.qname()) {
                    read_names.push(qname.to_string());
                }
            } else {
                break;
            }
        }

        // Verify correct reads are kept
        assert!(
            read_names.contains(&"read_keep".to_string()),
            "Read with lower subtraction score should be kept"
        );
        assert!(
            read_names.contains(&"read_unknown".to_string()),
            "Read not in subtraction should be kept"
        );

        // Verify eliminated reads are not present
        assert!(
            !read_names.contains(&"read_eliminate".to_string()),
            "Read with higher subtraction score should be eliminated"
        );
        assert!(
            !read_names.contains(&"read_equal".to_string()),
            "Read with equal scores should be eliminated"
        );

        // Verify the return value - should be 2 (read_eliminate and read_equal)
        assert_eq!(
            subtracted_count, 2,
            "Should have subtracted 2 reads (read_eliminate and read_equal)"
        );

        // Verify FASTQ output - should only contain read_keep and read_unknown
        let fastq_output = fs::read_to_string(&output_fastq_path).unwrap();
        assert!(
            fastq_output.contains("@read_keep"),
            "FASTQ output should contain read_keep"
        );
        assert!(
            fastq_output.contains("@read_unknown"),
            "FASTQ output should contain read_unknown"
        );
        assert!(
            !fastq_output.contains("@read_eliminate"),
            "FASTQ output should not contain read_eliminate"
        );
        assert!(
            !fastq_output.contains("@read_equal"),
            "FASTQ output should not contain read_equal"
        );
    }
}
