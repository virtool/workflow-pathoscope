use rust_htslib::{bam, bam::Read};
use std::collections::HashSet;
use std::io;
use std::path::Path;
use thiserror::Error;


#[derive(Error, Debug)]
pub enum StreamProcessorError {
    #[error("Failed to open SAM/BAM file '{path}': {source}")]
    FileOpen { path: String, source: io::Error },
    
    #[error("Failed to read BAM/SAM record: {source}")]
    BamRead { source: rust_htslib::errors::Error },
    
    #[error("Failed to parse alignment score from record")]
    AlignmentScoreParse,
    
    #[error("Invalid UTF-8 in reference ID")]
    InvalidRefId,
    
    #[error("Invalid UTF-8 in read ID")]
    InvalidReadId,
    
    #[error("Failed to execute bowtie2 command: {source}")]
    Bowtie2Execution { source: io::Error },
    
    #[error("Bowtie2 process failed with exit code {code}")]
    Bowtie2Failed { code: i32 },
    
}

/// Extract candidate OTU reference IDs from a SAM/BAM file
/// 
/// Processes alignments and collects reference IDs for reads that meet the score cutoff.
/// This replaces the Python line-by-line processing in map_default_isolates.
/// 
/// # Arguments
/// * `sam_path` - Path to the SAM/BAM file to process
/// * `p_score_cutoff` - Minimum score threshold (AS:i score + read length)
/// 
/// # Returns
/// HashSet of reference IDs that have reads meeting the score cutoff
pub fn extract_candidate_otus_from_sam_file<P: AsRef<Path>>(
    sam_path: P,
    p_score_cutoff: f64,
) -> Result<HashSet<String>, StreamProcessorError> {
    let path_str = sam_path.as_ref().to_string_lossy().to_string();
    
    let mut reader = bam::Reader::from_path(&sam_path)
        .map_err(|e| StreamProcessorError::FileOpen {
            path: path_str,
            source: io::Error::new(io::ErrorKind::Other, e),
        })?;

    extract_candidate_otus_from_reader(&mut reader, p_score_cutoff)
}

/// Extract candidate OTU reference IDs from a BAM/SAM reader
/// 
/// Core processing logic that works with any BAM reader (file, stream, etc.)
/// 
/// # Arguments
/// * `reader` - BAM reader to process records from
/// * `p_score_cutoff` - Minimum score threshold (AS:i score + read length)
/// 
/// # Returns
/// HashSet of reference IDs that have reads meeting the score cutoff
pub fn extract_candidate_otus_from_reader(
    reader: &mut bam::Reader,
    p_score_cutoff: f64,
) -> Result<HashSet<String>, StreamProcessorError> {
    let mut candidate_otus = HashSet::new();
    let header = reader.header().clone();

    for result in reader.records() {
        let record = result.map_err(|e| StreamProcessorError::BamRead { source: e })?;

        // Skip unmapped reads
        if record.is_unmapped() {
            continue;
        }

        // Get reference name - skip if unmapped (tid < 0)
        if record.tid() < 0 {
            continue;
        }

        let ref_name_bytes = header.tid2name(record.tid() as u32);
        let ref_id = std::str::from_utf8(ref_name_bytes)
            .map_err(|_| StreamProcessorError::InvalidRefId)?
            .to_string();

        // Skip if reference is "*" (unmapped)
        if ref_id == "*" {
            continue;
        }

        // Calculate total score (AS score + read length) - same logic as existing code
        let total_score = calculate_total_score(&record)?;

        // Apply score cutoff
        if total_score >= p_score_cutoff {
            candidate_otus.insert(ref_id);
        }
    }

    Ok(candidate_otus)
}

/// Calculate total score for a BAM record (AS score + read length)
/// 
/// Follows the same logic as the existing Python SamLine.score calculation
/// and the Rust parse_isolate_scores function.
/// 
/// # Arguments
/// * `record` - BAM record to calculate score for
/// 
/// # Returns
/// Total score (AS alignment score + read length)
fn calculate_total_score(record: &bam::Record) -> Result<f64, StreamProcessorError> {
    let read_length = record.seq_len() as f64;

    // Get alignment score from AS:i: auxiliary field
    let as_score = match record.aux(b"AS") {
        Ok(aux) => match aux {
            rust_htslib::bam::record::Aux::I32(score) => score as f64,
            rust_htslib::bam::record::Aux::I8(score) => score as f64,
            rust_htslib::bam::record::Aux::I16(score) => score as f64,
            rust_htslib::bam::record::Aux::U8(score) => score as f64,
            rust_htslib::bam::record::Aux::U16(score) => score as f64,
            rust_htslib::bam::record::Aux::U32(score) => score as f64,
            _ => return Err(StreamProcessorError::AlignmentScoreParse),
        },
        Err(_) => return Err(StreamProcessorError::AlignmentScoreParse),
    };

    Ok(as_score + read_length)
}

/// Extract candidate OTU reference IDs from SAM text data
/// 
/// This function parses SAM format data directly from bytes without using rust-htslib.
/// It provides a streaming alternative that doesn't require temporary files or unsafe code.
/// 
/// # Arguments
/// * `sam_bytes` - Raw SAM format data as bytes (typically from subprocess stdout)
/// * `p_score_cutoff` - Minimum score threshold (AS:i score + read length)
/// 
/// # Returns
/// HashSet of reference IDs that have reads meeting the score cutoff
pub fn extract_candidate_otus_from_bytes(
    sam_bytes: &[u8],
    p_score_cutoff: f64,
) -> Result<HashSet<String>, StreamProcessorError> {
    let mut candidate_otus = HashSet::new();
    let sam_text = std::str::from_utf8(sam_bytes)
        .map_err(|_| StreamProcessorError::InvalidRefId)?;

    for line in sam_text.lines() {
        // Skip header lines (start with @)
        if line.starts_with('@') {
            continue;
        }

        // Skip empty lines
        if line.trim().is_empty() {
            continue;
        }

        // Parse SAM line - tab-separated format
        let fields: Vec<&str> = line.split('\t').collect();
        
        // SAM format requires at least 11 fields
        if fields.len() < 11 {
            continue;
        }

        // Extract key fields:
        // 0: QNAME (read name)
        // 1: FLAG 
        // 2: RNAME (reference name)
        // 9: SEQ (read sequence)
        let flag: u16 = fields[1].parse().unwrap_or(4); // Default to unmapped if parse fails
        let ref_name = fields[2];
        let seq_len = fields[9].len() as f64;

        // Skip unmapped reads (flag & 4 != 0) or reads mapping to "*"
        if (flag & 4) != 0 || ref_name == "*" {
            continue;
        }

        // Find AS:i score in the optional fields (starting from field 11)
        let mut as_score: Option<f64> = None;
        for field in &fields[11..] {
            if field.starts_with("AS:i:") {
                if let Ok(score) = field[5..].parse::<i32>() {
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
                candidate_otus.insert(ref_name.to_string());
            }
        }
    }

    Ok(candidate_otus)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_candidate_otus_basic() {
        // Test with the basic SAM file
        let result = extract_candidate_otus_from_sam_file(
            "example/rust/test_basic.sam", 
            0.01
        ).unwrap();

        // Should find all 3 references since scores are:
        // read1: AS:i:45 + 50 = 95.0
        // read2: AS:i:25 + 30 = 55.0  
        // read3: AS:i:2 + 40 = 42.0
        assert_eq!(result.len(), 2, "Should find 2 unique references");
        assert!(result.contains("ref1"), "Should contain ref1");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidate_otus_with_cutoff() {
        // Test with higher cutoff to filter some results
        let result = extract_candidate_otus_from_sam_file(
            "example/rust/test_basic.sam", 
            50.0
        ).unwrap();

        // Only read1 (95.0) and read2 (55.0) should pass, read3 (42.0) should be filtered
        assert_eq!(result.len(), 2, "Should find 2 references with high cutoff");
        assert!(result.contains("ref1"), "Should contain ref1");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidate_otus_very_high_cutoff() {
        // Test with very high cutoff
        let result = extract_candidate_otus_from_sam_file(
            "example/rust/test_basic.sam", 
            100.0
        ).unwrap();

        // No reads should pass this cutoff
        assert_eq!(result.len(), 0, "Should find no references with very high cutoff");
    }

    #[test]
    fn test_extract_candidate_otus_file_not_found() {
        let result = extract_candidate_otus_from_sam_file(
            "/nonexistent/file.sam", 
            0.01
        );
        
        assert!(result.is_err(), "Should fail with non-existent file");
        match result.unwrap_err() {
            StreamProcessorError::FileOpen { path, .. } => {
                assert!(path.contains("nonexistent"), "Error should mention the file path");
            },
            _ => panic!("Expected FileOpen error"),
        }
    }

    #[test]
    fn test_calculate_total_score() {
        // This is an integration test that would require creating a mock BAM record
        // For now, we test the logic through the file processing tests above
        // The calculate_total_score function is tested indirectly
    }

    #[test]
    fn test_extract_candidate_otus_deduplication() {
        // Test that multiple reads mapping to the same reference only add it once
        let result = extract_candidate_otus_from_sam_file(
            "example/rust/test_basic.sam", 
            0.01
        ).unwrap();

        // Even if multiple reads map to ref1, it should only appear once in the set
        let ref1_count = result.iter().filter(|&r| r == "ref1").count();
        assert_eq!(ref1_count, 1, "Each reference should appear only once in the result set");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_basic() {
        // Test basic functionality with the same SAM data as test_basic.sam
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted
@SQ\tSN:ref1\tLN:1000
@SQ\tSN:ref2\tLN:2000
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45
read2\t0\tref2\t200\t255\t30M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t*\tAS:i:25
read3\t0\tref1\t300\t255\t40M\t*\t0\t0\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\t*\tAS:i:2";

        let result = extract_candidate_otus_from_bytes(sam_data, 0.01).unwrap();

        // Should find 2 unique references since scores are:
        // read1: AS:i:45 + 50 = 95.0
        // read2: AS:i:25 + 30 = 55.0  
        // read3: AS:i:2 + 40 = 42.0
        assert_eq!(result.len(), 2, "Should find 2 unique references");
        assert!(result.contains("ref1"), "Should contain ref1");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_with_cutoff() {
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted
@SQ\tSN:ref1\tLN:1000
@SQ\tSN:ref2\tLN:2000
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45
read2\t0\tref2\t200\t255\t30M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t*\tAS:i:25
read3\t0\tref1\t300\t255\t40M\t*\t0\t0\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\t*\tAS:i:2";

        let result = extract_candidate_otus_from_bytes(sam_data, 50.0).unwrap();

        // Only read1 (95.0) and read2 (55.0) should pass, read3 (42.0) should be filtered
        assert_eq!(result.len(), 2, "Should find 2 references with high cutoff");
        assert!(result.contains("ref1"), "Should contain ref1");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_very_high_cutoff() {
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted
@SQ\tSN:ref1\tLN:1000
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45";

        let result = extract_candidate_otus_from_bytes(sam_data, 100.0).unwrap();

        // No reads should pass this cutoff (read1 is 95.0)
        assert_eq!(result.len(), 0, "Should find no references with very high cutoff");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_unmapped_reads() {
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted
read1\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAA\t*
read2\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45";

        let result = extract_candidate_otus_from_bytes(sam_data, 0.01).unwrap();

        // Only mapped read should be counted
        assert_eq!(result.len(), 1, "Should ignore unmapped reads");
        assert!(result.contains("ref1"), "Should contain ref1");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_no_as_score() {
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*
read2\t0\tref2\t100\t255\t30M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t*\tAS:i:25";

        let result = extract_candidate_otus_from_bytes(sam_data, 0.01).unwrap();

        // Only read with AS score should be counted
        assert_eq!(result.len(), 1, "Should ignore reads without AS score");
        assert!(result.contains("ref2"), "Should contain ref2");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_empty_data() {
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted";

        let result = extract_candidate_otus_from_bytes(sam_data, 0.01).unwrap();

        assert_eq!(result.len(), 0, "Should handle empty SAM data");
    }

    #[test]
    fn test_extract_candidate_otus_from_bytes_malformed_lines() {
        let sam_data = b"@HD\tVN:1.0\tSO:unsorted
incomplete_line_with_few_fields\t0\tref1
read1\t0\tref1\t100\t255\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tAS:i:45";

        let result = extract_candidate_otus_from_bytes(sam_data, 0.01).unwrap();

        // Should skip malformed line and process valid one
        assert_eq!(result.len(), 1, "Should skip malformed lines");
        assert!(result.contains("ref1"), "Should contain ref1");
    }
}