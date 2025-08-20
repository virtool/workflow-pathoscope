use rust_htslib::{bam, bam::Read};
use std::{
    fmt::Debug,
    path::Path,
};

/// Supported alignment file formats
#[derive(Debug, Clone, PartialEq)]
pub enum AlignmentFormat {
    Sam,
    Bam,
}

/// Detect the format of an alignment file based on its extension
/// 
/// This function provides format detection for alignment files to ensure
/// proper handling by downstream consumers.
/// 
/// # Arguments
/// * `path` - Path to the alignment file
/// 
/// # Returns
/// The detected format (SAM or BAM)
pub fn detect_alignment_format<P: AsRef<Path>>(path: P) -> AlignmentFormat {
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("bam") => AlignmentFormat::Bam,
        Some("sam") => AlignmentFormat::Sam,
        _ => {
            // Default to SAM for unknown extensions or no extension
            // rust-htslib will handle format detection at the binary level
            AlignmentFormat::Sam
        }
    }
}

/// Minimal alignment data for memory-efficient storage
#[derive(Debug, Clone)]
pub struct MinimalAlignment {
    pub read_idx: i32,      // Index into reads vector
    pub ref_idx: i32,       // Index into refs vector  
    pub position: u32,      // Start position
    pub read_length: u16,   // Read length
}

/// Stores the desired fields of a .SAM record and the line itself as a String
#[derive(Debug, Clone)]
pub struct SamLine {
    pub read_id: String,
    pub read_length: usize,
    pub position: u32,
    pub score: Option<f64>,
    pub btws_flg: u32,
    pub unmapped: bool,
    pub ref_id: String,
    pub sam_fields: Vec<String>,
    pub line: String,
}



/// Parse an alignment file (SAM or BAM) from a file path and return all valid SamLine objects
/// 
/// This function uses rust-htslib which automatically detects and handles both SAM and BAM formats.
/// The format is detected at the binary level by rust-htslib, but we also provide extension-based
/// detection for clarity and debugging purposes.
/// 
/// # Arguments
/// * `alignment_path` - Path to the SAM or BAM file
/// * `p_score_cutoff` - Optional minimum score threshold for alignments
/// 
/// # Returns
/// Vector of SamLine objects that meet the score cutoff
pub fn parse_alignment<P: AsRef<Path>>(
    alignment_path: P,
    p_score_cutoff: Option<f64>,
) -> Result<Vec<SamLine>, String> {
    let p_score_cutoff = p_score_cutoff.unwrap_or(0.01);
    
    // Detect format for debugging/logging purposes
    let detected_format = detect_alignment_format(&alignment_path);
    
    let mut reader = bam::Reader::from_path(&alignment_path)
        .map_err(|e| format!("Failed to open alignment file '{}' (detected format: {:?}): {}", 
                            alignment_path.as_ref().display(), detected_format, e))?;

    let header = reader.header().clone();
    let mut sam_lines = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| format!("Error reading record: {}", e))?;

        // Skip unmapped reads
        if record.is_unmapped() {
            continue;
        }

        // Get read ID (qname)
        let read_id = std::str::from_utf8(record.qname())
            .map_err(|e| format!("Invalid UTF-8 in read ID: {}", e))?
            .to_string();

        // Get reference name
        let name_bytes = header.tid2name(record.tid() as u32);
        let ref_id = std::str::from_utf8(name_bytes).unwrap_or("*").to_string();

        // Get position (1-based in SAM format)
        let position = record.pos() as u32 + 1;

        // Get read length
        let read_length = record.seq_len();

        // Get alignment score from AS:i: auxiliary field
        let as_score = match record.aux(b"AS") {
            Ok(aux) => match aux {
                rust_htslib::bam::record::Aux::I32(score) => score as f64,
                rust_htslib::bam::record::Aux::I8(score) => score as f64,
                rust_htslib::bam::record::Aux::I16(score) => score as f64,
                rust_htslib::bam::record::Aux::U8(score) => score as f64,
                rust_htslib::bam::record::Aux::U16(score) => score as f64,
                rust_htslib::bam::record::Aux::U32(score) => score as f64,
                _ => continue, // Skip records with unexpected AS type
            },
            Err(_) => continue, // Skip records without AS field
        };

        // Calculate total score (AS score + read length, matching original logic)
        let total_score = as_score + read_length as f64;

        // Apply score cutoff
        if total_score > p_score_cutoff {
            let sam_line = SamLine {
                read_id,
                read_length,
                position,
                score: Some(total_score),
                btws_flg: record.flags() as u32,
                unmapped: record.is_unmapped(),
                ref_id,
                sam_fields: Vec::default(), // Not needed for htslib version
                line: String::default(),    // Not needed for htslib version
            };

            sam_lines.push(sam_line);
        }
    }

    Ok(sam_lines)
}

/// Legacy alias for parse_alignment to maintain backward compatibility
/// 
/// # Deprecated
/// Use `parse_alignment` instead as this function now supports both SAM and BAM formats
pub fn parse_sam<P: AsRef<Path>>(
    sam_path: P,
    p_score_cutoff: Option<f64>,
) -> Result<Vec<SamLine>, String> {
    parse_alignment(sam_path, p_score_cutoff)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_alignment_format() {
        assert_eq!(detect_alignment_format("test.sam"), AlignmentFormat::Sam);
        assert_eq!(detect_alignment_format("test.bam"), AlignmentFormat::Bam);
        assert_eq!(detect_alignment_format("test.txt"), AlignmentFormat::Sam); // defaults to SAM
        assert_eq!(detect_alignment_format("test"), AlignmentFormat::Sam); // no extension defaults to SAM
        
        // Test with Path objects
        use std::path::Path;
        assert_eq!(detect_alignment_format(Path::new("data/alignments.bam")), AlignmentFormat::Bam);
        assert_eq!(detect_alignment_format(Path::new("data/alignments.sam")), AlignmentFormat::Sam);
    }

    #[test]
    fn test_parse_sam_basic() {
        // Test with default cutoff (0.01) - should include all reads
        let sam_lines: Vec<SamLine> = parse_sam("example/rust/test_basic.sam", None)
            .unwrap();

        assert_eq!(sam_lines.len(), 3);

        // Check first read
        assert_eq!(sam_lines[0].read_id, "read1");
        assert_eq!(sam_lines[0].ref_id, "ref1");
        assert_eq!(sam_lines[0].position, 100);
        assert_eq!(sam_lines[0].read_length, 50);
        assert_eq!(sam_lines[0].score, Some(95.0)); // AS:i:45 + read_length:50

        // Check second read
        assert_eq!(sam_lines[1].read_id, "read2");
        assert_eq!(sam_lines[1].ref_id, "ref2");
        assert_eq!(sam_lines[1].position, 200);
        assert_eq!(sam_lines[1].read_length, 30);
        assert_eq!(sam_lines[1].score, Some(55.0)); // AS:i:25 + read_length:30

        // Check third read
        assert_eq!(sam_lines[2].read_id, "read3");
        assert_eq!(sam_lines[2].ref_id, "ref1");
        assert_eq!(sam_lines[2].position, 300);
        assert_eq!(sam_lines[2].read_length, 40);
        assert_eq!(sam_lines[2].score, Some(42.0)); // AS:i:2 + read_length:40
    }

    #[test]
    fn test_parse_sam_with_cutoff() {
        // Test with high cutoff (50.0) - should only include read1 and read2
        let sam_lines: Vec<SamLine> = parse_sam("example/rust/test_cutoff.sam", Some(50.0))
            .unwrap();

        assert_eq!(sam_lines.len(), 2);
        assert_eq!(sam_lines[0].read_id, "read1");
        assert_eq!(sam_lines[1].read_id, "read2");
    }

    #[test]
    fn test_parse_alignment_sam_format() {
        // Test that parse_alignment works with SAM files (same as parse_sam)
        let sam_lines: Vec<SamLine> = parse_alignment("example/rust/test_basic.sam", None)
            .unwrap();

        assert_eq!(sam_lines.len(), 3);
        assert_eq!(sam_lines[0].read_id, "read1");
        assert_eq!(sam_lines[0].ref_id, "ref1");
        assert_eq!(sam_lines[0].score, Some(95.0)); // AS:i:45 + read_length:50
    }

    #[test]
    fn test_parse_alignment_backward_compatibility() {
        // Ensure parse_sam and parse_alignment return identical results for SAM files
        let sam_results = parse_sam("example/rust/test_basic.sam", None).unwrap();
        let alignment_results = parse_alignment("example/rust/test_basic.sam", None).unwrap();
        
        assert_eq!(sam_results.len(), alignment_results.len());
        for (sam_line, align_line) in sam_results.iter().zip(alignment_results.iter()) {
            assert_eq!(sam_line.read_id, align_line.read_id);
            assert_eq!(sam_line.ref_id, align_line.ref_id);
            assert_eq!(sam_line.score, align_line.score);
            assert_eq!(sam_line.position, align_line.position);
        }
    }


    #[test]
    fn test_parse_sam_file_not_found() {
        let result = parse_sam("/nonexistent/file.sam", None);
        assert!(result.is_err());
        // Updated to match new error message format (includes "alignment file")
        assert!(result.unwrap_err().contains("Failed to open alignment file"));
    }
}