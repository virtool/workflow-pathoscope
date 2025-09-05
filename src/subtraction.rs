use crate::sam::{SamReader, extract_alignment_score};
use rust_htslib::bam;
use rust_htslib::bam::Format;
use rustc_hash::FxHashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use thiserror::Error;
use log::info;

#[derive(Error, Debug)]
pub enum BamProcessingError {
    #[error("Failed to create output file '{path}': {source}")]
    FileCreate { path: String, source: std::io::Error },
    
    #[error("Failed to parse SAM file: {0}")]
    SamParse(String),
    
    #[error("Failed to write BAM/SAM record: {source}")]
    BamWrite { source: rust_htslib::errors::Error },
    
    #[error("Failed to write to output file: {source}")]
    WriteOutput { source: std::io::Error },
}

#[derive(Debug, Clone)]
pub struct SubtractionProcessor {
    subtraction_scores: FxHashMap<String, f32>,
}


impl SubtractionProcessor {
    /// Creates a new SubtractionProcessor with the given subtraction scores
    pub fn new(subtraction_scores: FxHashMap<String, f32>) -> Self {
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

}

/// Parse subtraction BAM file and return a map of read IDs to scores.
pub fn parse_subtraction_scores_from_bam(path: &str) -> Result<FxHashMap<String, f32>, BamProcessingError> {
    info!("parsing subtraction SAM file: {}", path);
    
    let mut reader = SamReader::new(path)
        .map_err(BamProcessingError::SamParse)?;

    // Pre-allocate HashMap with estimated capacity
    let mut high_scores: FxHashMap<String, f32> = FxHashMap::with_capacity_and_hasher(100_000, Default::default());
    
    reader.stream_chunks(|chunk| {
        for record in chunk {
            // Skip unmapped reads
            if record.is_unmapped() {
                continue;
            }

            // Get read ID
            let read_id = std::str::from_utf8(record.qname())
                .map_err(|_| "Invalid UTF-8 in read ID")?
                .to_string();

            // Get alignment score
            if let Some(total_score) = extract_alignment_score(record) {
                high_scores.insert(read_id, total_score as f32);
            }
        }
        Ok(())
    }).map_err(BamProcessingError::SamParse)?;

    info!("parsed {} subtraction scores from {}", high_scores.len(), path);
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
) -> Result<(), BamProcessingError> {
    info!("starting subtraction elimination: isolate={}, subtraction={}, output={}",
          isolate_sam_path, subtraction_sam_path, output_sam_path);
    
    // Parse subtraction scores
    let subtraction_scores = parse_subtraction_scores_from_bam(subtraction_sam_path)?;
    let processor = SubtractionProcessor::new(subtraction_scores);
    
    // Process isolate file
    let subtracted_ids = process_isolate_file(isolate_sam_path, output_sam_path, &processor)?;
    
    info!("subtraction complete: {} reads eliminated", subtracted_ids.len());
    
    // Write subtracted IDs file
    write_subtracted_ids_file(output_sam_path, &subtracted_ids)?;
    
    Ok(())
}

/// Process the isolate SAM file and write filtered output using rust-htslib
pub fn process_isolate_file(
    input_path: &str,
    output_path: &str,
    processor: &SubtractionProcessor,
) -> Result<HashSet<String>, BamProcessingError> {
    
    let mut reader = SamReader::new(input_path)
        .map_err(|_| BamProcessingError::SamParse("Failed to open SAM file".to_string()))?;
    
    // Get the header from input to preserve SQ lines and other headers
    let header_view = reader.header().clone();
    
    // Create a new Header from the HeaderView
    let header = bam::Header::from_template(&header_view);
    
    // Create output writer with the same header
    let mut writer = bam::Writer::from_path(output_path, &header, Format::Bam)
        .map_err(|e| BamProcessingError::BamWrite { source: e })?;

    let mut subtracted_read_ids = HashSet::new();

    // Process each chunk
    reader.stream_chunks(|chunk| {
        for record in chunk {
            // Skip unmapped reads
            if record.is_unmapped() {
                continue;
            }
            
            // Get read ID as string slice first (no allocation)
            let read_id_str = std::str::from_utf8(record.qname())
                .map_err(|_| "Invalid UTF-8 in read ID")?;
            
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
            
            // Calculate alignment score using shared function
            let isolate_score = match extract_alignment_score(record) {
                Some(score) => score as f32,
                None => continue,
            };
            
            // Check if this read should be eliminated
            if processor.should_eliminate(read_id_str, isolate_score) {
                // Only allocate string when we need to store it
                subtracted_read_ids.insert(read_id_str.to_string());
            } else {
                // Write the record to output
                writer.write(record)
                    .map_err(|e| format!("Failed to write BAM record: {}", e))?;
            }
        }
        Ok(())
    }).map_err(|e| BamProcessingError::SamParse(e))?;

    Ok(subtracted_read_ids)
}

/// Write subtracted read IDs to file in same directory as output SAM
fn write_subtracted_ids_file(
    output_sam_path: &str,
    subtracted_ids: &HashSet<String>,
) -> Result<(), BamProcessingError> {
    // Create subtracted_read_ids.txt in the same directory as the output SAM file
    let output_path = Path::new(output_sam_path);
    let subtracted_ids_path = output_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("subtracted_read_ids.txt");

    let mut subtracted_read_ids_file = File::create(&subtracted_ids_path)
        .map_err(|e| BamProcessingError::FileCreate {
            path: subtracted_ids_path.to_string_lossy().to_string(),
            source: e,
        })?;

    for read_id in subtracted_ids {
        writeln!(&mut subtracted_read_ids_file, "{}", read_id)
            .map_err(|e| BamProcessingError::WriteOutput { source: e })?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_subtraction_scores() -> FxHashMap<String, f32> {
        let mut scores = FxHashMap::default();
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
    fn test_subtraction_processor_empty_scores() {
        let processor = SubtractionProcessor::new(FxHashMap::default());
        
        // With no subtraction scores, nothing should be eliminated
        assert!(!processor.should_eliminate("any_read", 1000.0));
        assert!(!processor.should_eliminate("any_read", 0.0));
    }


    #[test]
    fn test_parse_subtraction_sam() {
        let result = parse_subtraction_scores_from_bam("example/rust/test_basic.sam").unwrap();

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