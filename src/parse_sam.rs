use rust_htslib::{bam, bam::Read};
use std::{
    fmt::Debug,
    path::Path,
};

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



/// Parse a SAM file from a file path and return all valid SamLine objects
pub fn parse_sam<P: AsRef<Path>>(
    sam_path: P,
    p_score_cutoff: Option<f64>,
) -> Result<Vec<SamLine>, String> {
    let p_score_cutoff = p_score_cutoff.unwrap_or(0.01);
    
    let mut reader = bam::Reader::from_path(&sam_path)
        .map_err(|e| format!("Failed to open SAM file '{}': {}", sam_path.as_ref().display(), e))?;

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

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_parse_sam_file_not_found() {
        let result = parse_sam("/nonexistent/file.sam", None);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Failed to open SAM file"));
    }
}