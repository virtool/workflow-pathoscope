use crate::parse_sam::parse_sam;
use rust_htslib::{bam, bam::Read};
use rust_htslib::bam::Format;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SubtractionError {
    #[error("Failed to open SAM file '{path}': {source}")]
    FileOpen { path: String, source: std::io::Error },
    
    #[error("Failed to create output file '{path}': {source}")]
    FileCreate { path: String, source: std::io::Error },
    
    #[error("Failed to parse SAM file: {0}")]
    SamParse(String),
    
    #[error("Failed to read BAM/SAM record: {source}")]
    BamRead { source: rust_htslib::errors::Error },
    
    #[error("Failed to write BAM/SAM record: {source}")]
    BamWrite { source: rust_htslib::errors::Error },
    
    #[error("Insufficient fields in SAM line")]
    InsufficientFields,
    
    #[error("Failed to parse alignment score from field '{field}': {source}")]
    AlignmentScoreParse { field: String, source: std::num::ParseFloatError },
    
    #[error("Failed to read line from file: {source}")]
    LineRead { source: std::io::Error },
    
    #[error("Failed to write to output file: {source}")]
    WriteOutput { source: std::io::Error },
}

#[derive(Debug, Clone)]
pub struct SubtractionProcessor {
    subtraction_scores: HashMap<String, f32>,
}

#[derive(Debug, PartialEq)]
pub enum ProcessResult {
    Keep(String),       // SAM line to write to output (headers, comments, mapped reads)
    Eliminate(String),  // Read ID that was eliminated
}

impl SubtractionProcessor {
    /// Creates a new SubtractionProcessor with the given subtraction scores
    pub fn new(subtraction_scores: HashMap<String, f32>) -> Self {
        Self {
            subtraction_scores,
        }
    }

    /// Check if a read should be eliminated based on its isolate score vs subtraction score
    pub fn should_eliminate(&self, read_id: &str, isolate_score: f32) -> bool {
        match self.subtraction_scores.get(read_id) {
            Some(subtraction_score) => subtraction_score >= &isolate_score,
            None => false,
        }
    }

    /// Process a SAM line and determine what action to take
    /// Returns None for lines that should be skipped (empty, unmapped, malformed)
    pub fn process_sam_line(&self, line: &str) -> Result<Option<ProcessResult>, SubtractionError> {
        // Skip empty lines
        if line.trim().is_empty() {
            return Ok(None);
        }

        // Keep header lines (@) and comment lines (#)
        match line.chars().next() {
            Some('@') | Some('#') => return Ok(Some(ProcessResult::Keep(line.to_string()))),
            None => return Ok(None),
            _ => {}
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Skip lines with insufficient fields for a SAM record
        if fields.len() < 11 {
            return Ok(None);
        }

        let read_id = fields[0];
        let reference = fields[2];

        // Skip unmapped reads (reference is *)
        if reference == "*" {
            return Ok(None);
        }

        // Calculate alignment score
        let isolate_score = self.find_sam_align_score(&fields)?;

        // Check if this read should be eliminated
        if self.should_eliminate(read_id, isolate_score) {
            Ok(Some(ProcessResult::Eliminate(read_id.to_string())))
        } else {
            Ok(Some(ProcessResult::Keep(line.to_string())))
        }
    }

    /// Find the Pathoscope alignment score for a SAM line
    /// Returns AS score + read length
    fn find_sam_align_score(&self, fields: &[&str]) -> Result<f32, SubtractionError> {
        if fields.len() < 10 {
            return Err(SubtractionError::InsufficientFields);
        }

        let read_length = fields[9].chars().count() as f32;
        let mut a_score: f32 = 0.0;

        // Look for AS:i: tag in optional fields (starting from field 11)
        for field in fields.iter().skip(11) {
            if let Some(stripped) = field.strip_prefix("AS:i:") {
                a_score = stripped.parse().map_err(|e| {
                    SubtractionError::AlignmentScoreParse {
                        field: field.to_string(),
                        source: e,
                    }
                })?;
                break;
            }
        }

        Ok(a_score + read_length)
    }
}

/// Parse subtraction SAM file using parse_sam module and return scores for each read
pub fn parse_subtraction_sam(path: &str) -> Result<HashMap<String, f32>, SubtractionError> {
    let sam_lines = parse_sam(path, None)
        .map_err(SubtractionError::SamParse)?;

    let mut high_scores: HashMap<String, f32> = HashMap::new();
    
    for sam_line in sam_lines {
        if let Some(score) = sam_line.score {
            high_scores.insert(sam_line.read_id, score as f32);
        }
    }

    Ok(high_scores)
}

/// Main elimination function - pure Rust implementation
/// 
/// Processes isolate and subtraction SAM files to eliminate reads based on subtraction scores.
/// Creates both filtered SAM output and subtracted read IDs file.
pub fn eliminate_subtraction(
    isolate_sam_path: &str,
    subtraction_sam_path: &str,
    output_sam_path: &str,
) -> Result<(), SubtractionError> {
    // Parse subtraction scores
    let subtraction_scores = parse_subtraction_sam(subtraction_sam_path)?;
    let processor = SubtractionProcessor::new(subtraction_scores);
    
    // Process isolate file
    let subtracted_ids = process_isolate_file(isolate_sam_path, output_sam_path, &processor)?;
    
    // Write subtracted IDs file
    write_subtracted_ids_file(output_sam_path, &subtracted_ids)?;
    
    Ok(())
}

/// Process the isolate SAM file and write filtered output using rust-htslib
pub fn process_isolate_file(
    input_path: &str,
    output_path: &str,
    processor: &SubtractionProcessor,
) -> Result<HashSet<String>, SubtractionError> {
    // Open input file with rust-htslib reader
    let mut reader = bam::Reader::from_path(input_path)
        .map_err(|e| SubtractionError::BamRead { source: e })?;
    
    // Get the header from input to preserve SQ lines and other headers
    let header_view = reader.header().clone();
    
    // Create a new Header from the HeaderView
    let header = bam::Header::from_template(&header_view);
    
    // Create output writer with the same header
    let mut writer = bam::Writer::from_path(output_path, &header, Format::Bam)
        .map_err(|e| SubtractionError::BamWrite { source: e })?;

    let mut subtracted_read_ids = HashSet::new();

    // Process each record
    for result in reader.records() {
        let record = result.map_err(|e| SubtractionError::BamRead { source: e })?;
        
        // Skip unmapped reads
        if record.is_unmapped() {
            continue;
        }
        
        // Get read ID
        let read_id = std::str::from_utf8(record.qname())
            .map_err(|_| SubtractionError::InsufficientFields)?
            .to_string();
        
        // Get reference name
        let ref_name = if record.tid() >= 0 {
            std::str::from_utf8(header_view.tid2name(record.tid() as u32))
                .unwrap_or("*")
        } else {
            "*"
        };
        
        // Skip if reference is unmapped
        if ref_name == "*" {
            continue;
        }
        
        // Calculate alignment score (AS score + read length)
        let read_length = record.seq_len() as f32;
        let as_score = match record.aux(b"AS") {
            Ok(aux) => match aux {
                rust_htslib::bam::record::Aux::I32(score) => score as f32,
                rust_htslib::bam::record::Aux::I8(score) => score as f32,
                rust_htslib::bam::record::Aux::I16(score) => score as f32,
                rust_htslib::bam::record::Aux::U8(score) => score as f32,
                rust_htslib::bam::record::Aux::U16(score) => score as f32,
                rust_htslib::bam::record::Aux::U32(score) => score as f32,
                _ => 0.0,
            },
            Err(_) => 0.0,
        };
        let isolate_score = as_score + read_length;
        
        // Check if this read should be eliminated
        if processor.should_eliminate(&read_id, isolate_score) {
            subtracted_read_ids.insert(read_id);
        } else {
            // Write the record to output
            writer.write(&record)
                .map_err(|e| SubtractionError::BamWrite { source: e })?;
        }
    }

    Ok(subtracted_read_ids)
}

/// Write subtracted read IDs to file in same directory as output SAM
fn write_subtracted_ids_file(
    output_sam_path: &str,
    subtracted_ids: &HashSet<String>,
) -> Result<(), SubtractionError> {
    // Create subtracted_read_ids.txt in the same directory as the output SAM file
    let output_path = Path::new(output_sam_path);
    let subtracted_ids_path = output_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("subtracted_read_ids.txt");

    let mut subtracted_read_ids_file = File::create(&subtracted_ids_path)
        .map_err(|e| SubtractionError::FileCreate {
            path: subtracted_ids_path.to_string_lossy().to_string(),
            source: e,
        })?;

    for read_id in subtracted_ids {
        writeln!(&mut subtracted_read_ids_file, "{}", read_id)
            .map_err(|e| SubtractionError::WriteOutput { source: e })?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_subtraction_scores() -> HashMap<String, f32> {
        let mut scores = HashMap::new();
        scores.insert("read1".to_string(), 250.0);
        scores.insert("read2".to_string(), 300.0);
        scores.insert("read3".to_string(), 150.0);
        scores
    }

    #[test]
    fn test_should_eliminate_higher_subtraction_score() {
        let processor = SubtractionProcessor::new(create_test_subtraction_scores());
        
        // Isolate score (200) < subtraction score (250) -> eliminate
        assert!(processor.should_eliminate("read1", 200.0));
    }

    #[test]
    fn test_should_eliminate_equal_scores() {
        let processor = SubtractionProcessor::new(create_test_subtraction_scores());
        
        // Isolate score (250) == subtraction score (250) -> eliminate
        assert!(processor.should_eliminate("read1", 250.0));
    }

    #[test]
    fn test_should_not_eliminate_lower_subtraction_score() {
        let processor = SubtractionProcessor::new(create_test_subtraction_scores());
        
        // Isolate score (350) > subtraction score (300) -> keep
        assert!(!processor.should_eliminate("read2", 350.0));
    }

    #[test]
    fn test_should_not_eliminate_unknown_read() {
        let processor = SubtractionProcessor::new(create_test_subtraction_scores());
        
        // Read not in subtraction -> keep
        assert!(!processor.should_eliminate("unknown_read", 100.0));
    }

    #[test]
    fn test_process_sam_line_header() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let line = "@HD\tVN:1.6\tSO:unsorted";
        
        assert_eq!(
            processor.process_sam_line(line).unwrap(),
            Some(ProcessResult::Keep(line.to_string()))
        );
    }

    #[test]
    fn test_process_sam_line_comment() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let line = "# This is a comment";
        
        assert_eq!(
            processor.process_sam_line(line).unwrap(),
            Some(ProcessResult::Keep(line.to_string()))
        );
    }

    #[test]
    fn test_process_sam_line_empty() {
        let processor = SubtractionProcessor::new(HashMap::new());
        
        assert_eq!(processor.process_sam_line("").unwrap(), None);
        assert_eq!(processor.process_sam_line("   ").unwrap(), None);
    }

    #[test]
    fn test_process_sam_line_insufficient_fields() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let line = "read1\t0\tref1";
        
        assert_eq!(processor.process_sam_line(line).unwrap(), None);
    }

    #[test]
    fn test_process_sam_line_unmapped() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let line = "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tAS:i:50";
        
        assert_eq!(processor.process_sam_line(line).unwrap(), None);
    }

    #[test]
    fn test_process_sam_line_keep() {
        let mut scores = HashMap::new();
        scores.insert("read1".to_string(), 100.0);
        let processor = SubtractionProcessor::new(scores);
        
        // AS:i:150 + read_length:4 = 154.0 > subtraction_score:100.0 -> keep
        let line = "read1\t0\tref1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:150";
        
        match processor.process_sam_line(line).unwrap() {
            Some(ProcessResult::Keep(kept_line)) => assert_eq!(kept_line, line),
            other => panic!("Expected Some(Keep), got {:?}", other),
        }
    }

    #[test]
    fn test_process_sam_line_eliminate() {
        let mut scores = HashMap::new();
        scores.insert("read1".to_string(), 200.0);
        let processor = SubtractionProcessor::new(scores);
        
        // AS:i:150 + read_length:4 = 154.0 < subtraction_score:200.0 -> eliminate
        let line = "read1\t0\tref1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:150";
        
        match processor.process_sam_line(line).unwrap() {
            Some(ProcessResult::Eliminate(read_id)) => assert_eq!(read_id, "read1"),
            other => panic!("Expected Some(Eliminate), got {:?}", other),
        }
    }

    #[test]
    fn test_find_sam_align_score_basic() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let fields = vec!["read1", "0", "ref1", "100", "60", "4M", "*", "0", "0", "ACGT", "IIII", "AS:i:150"];
        
        let score = processor.find_sam_align_score(&fields).unwrap();
        assert_eq!(score, 154.0); // 150 + 4
    }

    #[test]
    fn test_find_sam_align_score_no_as_tag() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let fields = vec!["read1", "0", "ref1", "100", "60", "4M", "*", "0", "0", "ACGT", "IIII"];
        
        let score = processor.find_sam_align_score(&fields).unwrap();
        assert_eq!(score, 4.0); // 0 + 4 (no AS tag found)
    }

    #[test]
    fn test_find_sam_align_score_multiple_tags() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let fields = vec!["read1", "0", "ref1", "100", "60", "4M", "*", "0", "0", "ACGT", "IIII", "XM:i:0", "AS:i:150", "NM:i:0"];
        
        let score = processor.find_sam_align_score(&fields).unwrap();
        assert_eq!(score, 154.0); // 150 + 4
    }

    #[test]
    fn test_find_sam_align_score_insufficient_fields() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let fields = vec!["read1", "0", "ref1"];
        
        assert!(processor.find_sam_align_score(&fields).is_err());
    }

    #[test]
    fn test_find_sam_align_score_invalid_as_value() {
        let processor = SubtractionProcessor::new(HashMap::new());
        let fields = vec!["read1", "0", "ref1", "100", "60", "4M", "*", "0", "0", "ACGT", "IIII", "AS:i:invalid"];
        
        assert!(processor.find_sam_align_score(&fields).is_err());
    }

    #[test]
    fn test_process_sam_line_unknown_read() {
        let processor = SubtractionProcessor::new(create_test_subtraction_scores());
        let line = "unknown_read\t0\tref1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:50";
        
        // Unknown read should be kept regardless of score
        match processor.process_sam_line(line).unwrap() {
            Some(ProcessResult::Keep(kept_line)) => assert_eq!(kept_line, line),
            other => panic!("Expected Some(Keep), got {:?}", other),
        }
    }

    #[test]
    fn test_subtraction_processor_empty_scores() {
        let processor = SubtractionProcessor::new(HashMap::new());
        
        // With no subtraction scores, nothing should be eliminated
        assert!(!processor.should_eliminate("any_read", 1000.0));
        assert!(!processor.should_eliminate("any_read", 0.0));
    }

    #[test]
    fn test_process_sam_line_long_read_sequence() {
        let mut scores = HashMap::new();
        scores.insert("long_read".to_string(), 200.0);
        let processor = SubtractionProcessor::new(scores);
        
        let long_seq = "A".repeat(150);
        let long_qual = "I".repeat(150);
        let line = format!("long_read\t0\tref1\t100\t60\t150M\t*\t0\t0\t{}\t{}\tAS:i:50", long_seq, long_qual);
        
        // AS:i:50 + read_length:150 = 200.0 == subtraction_score:200.0 -> eliminate
        match processor.process_sam_line(&line).unwrap() {
            Some(ProcessResult::Eliminate(read_id)) => assert_eq!(read_id, "long_read"),
            other => panic!("Expected Some(Eliminate), got {:?}", other),
        }
    }

    #[test]
    fn test_parse_subtraction_sam() {
        let result = parse_subtraction_sam("example/rust/test_basic.sam").unwrap();

        // Expected scores based on the SAM file content:
        // read1: AS:i:45 + read_length:50 = 95.0
        // read2: AS:i:25 + read_length:30 = 55.0  
        // read3: AS:i:2 + read_length:40 = 42.0

        assert_eq!(result.len(), 3, "Should parse exactly 3 valid alignments");

        assert_eq!(
            result.get("read1"),
            Some(&95.0),
            "First read should have score 95.0"
        );

        assert_eq!(
            result.get("read2"),
            Some(&55.0),
            "Second read should have score 55.0"
        );

        assert_eq!(
            result.get("read3"),
            Some(&42.0),
            "Third read should have score 42.0"
        );
    }
}