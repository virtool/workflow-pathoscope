use std::collections::HashSet;
use log::info;
use pyo3::prelude::*;
use pyo3::exceptions::PyIOError;




/// Parse a single SAM line and extract candidate OTU information
/// 
/// This function processes one SAM line and determines if the read meets the score cutoff.
/// Used for testing and by the streaming functions.
/// 
/// # Arguments
/// * `line` - A SAM format line as string
/// * `p_score_cutoff` - Minimum score threshold (AS:i score + read length)
/// 
/// # Returns
/// Option containing the reference name if the read meets the cutoff, None otherwise
pub fn parse_sam_line(line: &str, p_score_cutoff: f64) -> Option<String> {
    // Skip header lines and empty lines
    if line.starts_with('@') || line.trim().is_empty() {
        return None;
    }

    // Parse SAM line - tab-separated format
    let fields: Vec<&str> = line.split('\t').collect();
    
    // SAM format requires at least 11 fields
    if fields.len() < 11 {
        return None;
    }

    // Extract key fields:
    // 1: FLAG 
    // 2: RNAME (reference name)
    // 9: SEQ (read sequence)
    let flag: u16 = fields[1].parse().unwrap_or(4); // Default to unmapped if parse fails
    let ref_name = fields[2];
    let seq_len = fields[9].len() as f64;

    // Skip unmapped reads (flag & 4 != 0) or reads mapping to "*"
    if (flag & 4) != 0 || ref_name == "*" {
        return None;
    }

    // Find AS:i score in the optional fields (starting from field 11)
    let mut as_score: Option<f64> = None;
    for field in &fields[11..] {
        if let Some(stripped) = field.strip_prefix("AS:i:") {
            if let Ok(score) = stripped.parse::<i32>() {
                as_score = Some(score as f64);
                break;
            }
        }
    }

    // Skip if no AS score found
    if let Some(as_score) = as_score {
        // Calculate total score (AS score + read length)
        let total_score = as_score + seq_len;

        // Apply score cutoff
        if total_score >= p_score_cutoff {
            return Some(ref_name.to_string());
        }
    }

    None
}

/// Extract candidate OTU reference IDs from SAM text data
/// 
/// This function parses SAM format data directly from text without using rust-htslib.
/// Used for testing and can be called by other functions that need SAM text parsing.
/// 
/// # Arguments
/// * `sam_text` - Raw SAM format data as string
/// * `p_score_cutoff` - Minimum score threshold (AS:i score + read length)
/// 
/// # Returns
/// HashSet of reference IDs that have reads meeting the score cutoff
pub fn extract_candidates_from_sam_text(
    sam_text: &str,
    p_score_cutoff: f64,
) -> HashSet<String> {
    let mut candidate_otus = HashSet::new();

    for line in sam_text.lines() {
        if let Some(ref_name) = parse_sam_line(line, p_score_cutoff) {
            candidate_otus.insert(ref_name);
        }
    }

    candidate_otus
}

/// Extract candidate OTU reference IDs by running bowtie2 directly with streaming
/// 
/// This function spawns a bowtie2 process directly from Rust and streams its output
/// to avoid memory issues with large SAM files. It processes SAM lines as they arrive
/// and returns only the unique reference IDs that meet the score cutoff.
/// 
/// # Arguments
/// * `bowtie_index_path` - Path to the bowtie2 index
/// * `read_paths` - List of paths to the input read files
/// * `proc` - Number of processor threads for bowtie2
/// * `p_score_cutoff` - Minimum score threshold (AS:i score + read length)
/// 
/// # Returns
/// Set of reference IDs that have reads meeting the score cutoff
pub fn find_candidate_otus_with_bowtie2(
    py: Python,
    bowtie_index_path: &str,
    read_paths: Vec<String>,
    proc: i32,
    p_score_cutoff: f64,
) -> PyResult<HashSet<String>> {
    use std::process::{Command, Stdio};
    use std::io::{BufRead, BufReader};
    
    info!("running bowtie2 directly from rust: index={}, reads={:?}, cutoff={}", 
          bowtie_index_path, read_paths, p_score_cutoff);
    py.allow_threads(|| {
        let mut cmd = Command::new("bowtie2");
        cmd.arg("-p").arg(proc.to_string())
           .arg("--local")
           .arg("--no-unal")
           .arg("--score-min").arg("L,20,1.0")
           .arg("-N").arg("0")
           .arg("-L").arg("15")
           .arg("-x").arg(bowtie_index_path)
           .arg("-U").arg(read_paths.join(","))
           .stdout(Stdio::piped())
           .stderr(Stdio::piped());
           
        info!("spawning bowtie2 process");
        let mut child = cmd.spawn()
            .map_err(|e| PyErr::new::<PyIOError, _>(format!("Failed to spawn bowtie2: {}", e)))?;
        
        let stdout = child.stdout.take().unwrap();
        let reader = BufReader::new(stdout);
        
        let mut candidate_otus = HashSet::new();
        let mut line_count = 0u64;
        let mut passing_count = 0u64;
        
        for line_result in reader.lines() {
            let line = line_result
                .map_err(|e| PyErr::new::<PyIOError, _>(format!("Error reading bowtie2 output: {}", e)))?;
            
            line_count += 1;
            
            // Use the extracted SAM parsing function
            if let Some(ref_name) = parse_sam_line(&line, p_score_cutoff) {
                candidate_otus.insert(ref_name);
                passing_count += 1;
            }
        }
        
        // Wait for bowtie2 to finish and check exit status
        let status = child.wait()
            .map_err(|e| PyErr::new::<PyIOError, _>(format!("Error waiting for bowtie2: {}", e)))?;
            
        if !status.success() {
            // Read stderr for error details
            let stderr_output = if let Some(mut stderr) = child.stderr.take() {
                let mut buf = String::new();
                let _ = std::io::Read::read_to_string(&mut stderr, &mut buf);
                buf
            } else {
                "Unknown error".to_string()
            };
            
            return Err(PyErr::new::<PyIOError, _>(format!(
                "bowtie2 failed with exit code {:?}: {}", 
                status.code(), 
                stderr_output
            )));
        }
        
        info!("processed {} sam lines, {} passed cutoff, found {} unique otus", 
              line_count, passing_count, candidate_otus.len());
        
        Ok(candidate_otus)
    })
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sam_line_basic() {
        let line = "read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45";
        let result = parse_sam_line(line, 0.01);
        
        // AS:i:45 + seq_len(50) = 95.0, should pass cutoff of 0.01
        assert_eq!(result, Some("ref1".to_string()));
    }

    #[test]
    fn test_parse_sam_line_below_cutoff() {
        let line = "read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45";
        let result = parse_sam_line(line, 100.0);
        
        // AS:i:45 + seq_len(50) = 95.0, should not pass cutoff of 100.0
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_sam_line_unmapped() {
        let line = "read1\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAA\t*";
        let result = parse_sam_line(line, 0.01);
        
        // Unmapped read (flag & 4 != 0), should return None
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_sam_line_no_as_score() {
        let line = "read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*";
        let result = parse_sam_line(line, 0.01);
        
        // No AS:i score, should return None
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_sam_line_header() {
        let line = "@HD\tVN:1.0\tSO:unsorted";
        let result = parse_sam_line(line, 0.01);
        
        // Header line, should return None
        assert_eq!(result, None);
    }

    #[test]
    fn test_extract_candidates_from_sam_text_basic() {
        let sam_data = "@HD\tVN:1.0\tSO:unsorted
@SQ\tSN:ref1\tLN:1000
@SQ\tSN:ref2\tLN:2000
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45
read2\t0\tref2\t200\t255\t30M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t*\tAS:i:25
read3\t0\tref1\t300\t255\t40M\t*\t0\t0\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\t*\tAS:i:2";

        let result = extract_candidates_from_sam_text(sam_data, 0.01);

        // Should find 2 unique references since scores are:
        // read1: AS:i:45 + 50 = 95.0
        // read2: AS:i:25 + 30 = 55.0  
        // read3: AS:i:2 + 40 = 42.0
        assert_eq!(result.len(), 2, "Should find 2 unique references");
        assert!(result.contains("ref1"), "Should contain ref1");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidates_from_sam_text_with_cutoff() {
        let sam_data = "@HD\tVN:1.0\tSO:unsorted
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45
read2\t0\tref2\t200\t255\t30M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t*\tAS:i:25
read3\t0\tref1\t300\t255\t40M\t*\t0\t0\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\t*\tAS:i:2";

        let result = extract_candidates_from_sam_text(sam_data, 50.0);

        // Only read1 (95.0) and read2 (55.0) should pass, read3 (42.0) should be filtered
        assert_eq!(result.len(), 2, "Should find 2 references with high cutoff");
        assert!(result.contains("ref1"), "Should contain ref1");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidates_from_sam_text_very_high_cutoff() {
        let sam_data = "@HD\tVN:1.0\tSO:unsorted
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45";

        let result = extract_candidates_from_sam_text(sam_data, 100.0);

        // No reads should pass this cutoff (read1 is 95.0)
        assert_eq!(result.len(), 0, "Should find no references with very high cutoff");
    }

    #[test]
    fn test_extract_candidates_from_sam_text_deduplication() {
        let sam_data = "@HD\tVN:1.0\tSO:unsorted
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45
read2\t0\tref1\t200\t255\t30M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t*\tAS:i:25";

        let result = extract_candidates_from_sam_text(sam_data, 0.01);

        // Even if multiple reads map to ref1, it should only appear once in the set
        assert_eq!(result.len(), 1, "Should deduplicate reference names");
        assert!(result.contains("ref1"), "Should contain ref1");
        let ref1_count = result.iter().filter(|&r| r == "ref1").count();
        assert_eq!(ref1_count, 1, "Each reference should appear only once in the result set");
    }

}