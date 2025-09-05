use crate::sam::{MinimalAlignment, SamReader, extract_alignment_score};
use crate::{UniqueReads, MultiMappingReads, MatrixResult};
use rustc_hash::FxHashMap;
use rust_htslib::{bam, bam::HeaderView};
use log::info;

/// A matrix containing alignment data and metadata.
/// 
/// This struct encapsulates all the data produced by matrix building,
/// providing a clean interface and methods for processing alignment data.
#[derive(Debug, Clone)]
pub struct PathoscopeMatrix {
    pub unique_reads: UniqueReads,
    pub multi_mapping_reads: MultiMappingReads,
    pub refs: Vec<String>,
    pub reads: Vec<String>,
    pub alignments: Vec<MinimalAlignment>,
    pub max_score: f64,
    pub min_score: f64,
}

impl PathoscopeMatrix {
    /// Create PathoscopeMatrix from raw alignment data
    /// 
    /// # Arguments
    /// * `read_alignments` - Map of read indices to their alignment data
    /// * `refs` - Reference sequence names
    /// * `reads` - Read sequence names  
    /// * `alignments` - Minimal alignment data for coverage calculation
    /// * `max_score` - Maximum alignment score found
    /// * `min_score` - Minimum alignment score found
    pub fn from_alignments(
        read_alignments: FxHashMap<i32, Vec<(i32, f64)>>,
        refs: Vec<String>,
        reads: Vec<String>,
        alignments: Vec<MinimalAlignment>,
        max_score: f64,
        min_score: f64,
    ) -> Self {
        // Classify reads as unique or multi-mapping
        let mut u_temp: MultiMappingReads = FxHashMap::default();
        let mut nu: MultiMappingReads = FxHashMap::default();

        for (read_index, read_alignments) in read_alignments {
            if read_alignments.len() == 1 {
                // Unique read: maps to exactly one reference
                let (ref_index, score) = read_alignments[0];
                u_temp.insert(
                    read_index,
                    (vec![ref_index], vec![score], vec![score], score),
                );
            } else {
                // Non-unique read: maps to multiple references
                let ref_indices: Vec<i32> = read_alignments.iter().map(|(ref_idx, _)| *ref_idx).collect();
                let scores: Vec<f64> = read_alignments.iter().map(|(_, score)| *score).collect();
                let max_score = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                
                nu.insert(
                    read_index,
                    (ref_indices, scores, vec![], max_score),
                );
            }
        }

        let mut matrix = PathoscopeMatrix {
            unique_reads: FxHashMap::default(),
            multi_mapping_reads: nu,
            refs,
            reads,
            alignments,
            max_score,
            min_score,
        };

        // Rescale scores and build final data structures
        matrix.rescale_scores(u_temp);
        matrix.build_unique_map();
        matrix.normalize_multi_mapping();

        let unique_count = matrix.unique_reads.len();
        let multi_count = matrix.multi_mapping_reads.len();
        info!("matrix created: {} unique reads, {} multi-mapping reads, score range [{:.2}, {:.2}]",
              unique_count, multi_count, min_score, max_score);

        matrix
    }

    /// Rescale alignment scores using the existing rescale_samscore function
    fn rescale_scores(&mut self, u_temp: MultiMappingReads) {
        let (u, nu) = rescale_samscore(u_temp, self.multi_mapping_reads.clone(), self.max_score, self.min_score);
        self.multi_mapping_reads = nu;
        
        // Convert u to unique_reads format and store for build_unique_map
        self.unique_reads = FxHashMap::default();
        for (read_idx, (ref_indices, scores, _, _)) in u {
            if let (Some(first_ref), Some(first_score)) = (ref_indices.first(), scores.first()) {
                self.unique_reads.insert(read_idx, (*first_ref, *first_score));
            }
        }
    }

    /// Build the final unique reads map (placeholder as rescale_scores now handles this)
    fn build_unique_map(&mut self) {
        // This is now handled in rescale_scores method
    }

    /// Normalize multi-mapping read scores so they sum to 1.0
    fn normalize_multi_mapping(&mut self) {
        let nu_keys: Vec<i32> = self.multi_mapping_reads.keys().cloned().collect();
        for k in nu_keys {
            if let Some(nu_entry) = self.multi_mapping_reads.get(&k) {
                let p_score_sum = nu_entry.1.iter().sum::<f64>();
                if let Some(nu_entry_mut) = self.multi_mapping_reads.get_mut(&k) {
                    nu_entry_mut.2 = nu_entry_mut
                        .1
                        .iter()
                        .map(|data| data / p_score_sum)
                        .collect();
                }
            }
        }
    }

    /// Convert to the legacy MatrixResult tuple format for backward compatibility
    pub fn into_matrix_result(self) -> MatrixResult {
        (
            self.unique_reads,
            self.multi_mapping_reads,
            self.refs,
            self.reads,
            self.alignments,
        )
    }
}

/// Process a single BAM record and extract alignment data
/// 
/// # Arguments
/// * `record` - BAM record to process
/// * `header` - BAM file header for reference name lookup
/// * `p_score_cutoff` - Minimum score threshold
/// * `h_read_id` - Mutable reference to read ID to index mapping
/// * `h_ref_id` - Mutable reference to reference ID to index mapping
/// * `refs` - Mutable reference to reference names vector
/// * `reads` - Mutable reference to read names vector
/// * `ref_count` - Mutable reference to reference counter
/// * `read_count` - Mutable reference to read counter
/// * `max_score` - Mutable reference to maximum score tracker
/// * `min_score` - Mutable reference to minimum score tracker
/// 
/// # Returns
/// Option containing (read_index, ref_index, minimal_alignment, score) if record is valid
fn process_bam_record(
    record: &bam::Record,
    header: &HeaderView,
    p_score_cutoff: f64,
    h_read_id: &mut FxHashMap<String, i32>,
    h_ref_id: &mut FxHashMap<String, i32>,
    refs: &mut Vec<String>,
    reads: &mut Vec<String>,
    ref_count: &mut i32,
    read_count: &mut i32,
    max_score: &mut f64,
    min_score: &mut f64,
) -> Option<(i32, i32, MinimalAlignment, f64)> {
    // Skip unmapped reads
    if record.is_unmapped() {
        return None;
    }

    // Get read ID (qname)
    let read_id = std::str::from_utf8(record.qname()).ok()?.to_string();

    // Get reference name
    let name_bytes = header.tid2name(record.tid() as u32);
    let ref_id = std::str::from_utf8(name_bytes).unwrap_or("*").to_string();

    // Get position (1-based in SAM format)
    let position = record.pos() as u32 + 1;

    // Get read length
    let read_length = record.seq_len();

    // Get alignment score using shared function
    let total_score = extract_alignment_score(record)?;

    // Apply score cutoff
    if total_score <= p_score_cutoff {
        return None;
    }

    // Track score range
    *min_score = total_score.min(*min_score);
    *max_score = total_score.max(*max_score);

    // Get or create reference index
    let ref_index = *h_ref_id.get(&ref_id).unwrap_or(&-1);
    let ref_index = if ref_index == -1 {
        let new_ref_index = *ref_count;
        h_ref_id.insert(ref_id.clone(), new_ref_index);
        refs.push(ref_id);
        *ref_count += 1;
        new_ref_index
    } else {
        ref_index
    };

    // Get or create read index
    let read_index = *h_read_id.get(&read_id).unwrap_or(&-1);
    let read_index = if read_index == -1 {
        let new_read_index = *read_count;
        h_read_id.insert(read_id.clone(), new_read_index);
        reads.push(read_id);
        *read_count += 1;
        new_read_index
    } else {
        read_index
    };

    // Create minimal alignment data
    let minimal_alignment = MinimalAlignment {
        read_idx: read_index,
        ref_idx: ref_index,
        position,
        read_length: read_length as u16,
    };

    Some((read_index, ref_index, minimal_alignment, total_score))
}

/// Build read alignments map by processing BAM file in chunks
/// 
/// # Arguments
/// * `reader` - SamReader for streaming
/// * `header` - BAM file header  
/// * `p_score_cutoff` - Minimum score threshold
/// 
/// # Returns
/// Tuple containing (read_alignments, refs, reads, minimal_alignments, max_score, min_score)
fn create_read_alignments_map(
    mut reader: SamReader,
    header: &HeaderView,
    p_score_cutoff: f64,
) -> Result<(FxHashMap<i32, Vec<(i32, f64)>>, Vec<String>, Vec<String>, Vec<MinimalAlignment>, f64, f64), String> {
    let mut h_read_id: FxHashMap<String, i32> = FxHashMap::default();
    let mut h_ref_id: FxHashMap<String, i32> = FxHashMap::default();
    let mut refs: Vec<String> = Vec::new();
    let mut reads: Vec<String> = Vec::new();
    let mut ref_count: i32 = 0;
    let mut read_count: i32 = 0;
    let mut max_score: f64 = 0.0;
    let mut min_score: f64 = 0.0;
    let mut minimal_alignments = Vec::new();
    let mut read_alignments: FxHashMap<i32, Vec<(i32, f64)>> = FxHashMap::default();

    reader.stream_chunks(|chunk| {
        // Process this chunk
        for record in chunk {
            if let Some((read_index, ref_index, minimal_alignment, total_score)) = process_bam_record(
                record,
                header,
                p_score_cutoff,
                &mut h_read_id,
                &mut h_ref_id,
                &mut refs,
                &mut reads,
                &mut ref_count,
                &mut read_count,
                &mut max_score,
                &mut min_score,
            ) {
                minimal_alignments.push(minimal_alignment);

                // Add alignment to read's collection (skip duplicates)
                let alignments = read_alignments.entry(read_index).or_default();
                if !alignments.iter().any(|(ref_idx, _)| *ref_idx == ref_index) {
                    alignments.push((ref_index, total_score));
                }
            }
        }
        Ok(())
    })?;

    Ok((read_alignments, refs, reads, minimal_alignments, max_score, min_score))
}

/// Build the EM matrix.
/// 
/// # Arguments
/// * `alignment_path` - Path to the SAM/BAM file
/// * `p_score_cutoff` - Optional score cutoff for alignments
pub fn build_matrix(
    alignment_path: &str,
    p_score_cutoff: Option<f64>,
) -> Result<MatrixResult, String> {
    let p_score_cutoff = p_score_cutoff.unwrap_or(0.01);
    
    info!("building matrix from '{}' with score cutoff {}",
          alignment_path, p_score_cutoff);
    
    let reader = SamReader::new(alignment_path)?;
    let header = reader.header().clone();
    
    // Build read alignments map using helper function
    let (read_alignments, refs, reads, minimal_alignments, max_score, min_score) = 
        create_read_alignments_map(reader, &header, p_score_cutoff)?;

    // Create PathoscopeMatrix and convert to legacy format
    let matrix = PathoscopeMatrix::from_alignments(
        read_alignments,
        refs,
        reads,
        minimal_alignments,
        max_score,
        min_score,
    );

    Ok(matrix.into_matrix_result())
}

/// modifies the scores of u and nu with respect to max_score and min_score
fn rescale_samscore(
    mut u: MultiMappingReads,
    mut nu: MultiMappingReads,
    max_score: f64,
    min_score: f64,
) -> (
    MultiMappingReads,
    MultiMappingReads,
) {
    let scaling_factor: f64 = if min_score < 0.0 {
        100.0 / (max_score - min_score)
    } else {
        100.0 / max_score
    };

    let u_keys: Vec<i32> = u.keys().cloned().collect();
    for k in u_keys {
        if let Some(entry) = u.get_mut(&k) {
            if min_score < 0.0 {
                entry.1[0] -= min_score;
            }
            entry.1[0] = f64::exp(entry.1[0] * scaling_factor);
            entry.3 = entry.1[0];
        }
    }

    let nu_keys: Vec<i32> = nu.keys().cloned().collect();
    for k in nu_keys {
        if let Some(entry) = nu.get_mut(&k) {
            entry.3 = 0.0;

            for i in 0..entry.1.len() {
                if min_score < 0.0 {
                    entry.1[i] -= min_score;
                }

                entry.1[i] = f64::exp(entry.1[i] * scaling_factor);

                if entry.1[i] > entry.3 {
                    entry.3 = entry.1[i];
                }
            }
        }
    }
    (u, nu)
}

#[cfg(test)]
mod tests {
    #![allow(unused)]

    use crate::matrix::*;
    use crate::*;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    use std::io::Read;

    extern crate yaml_rust;
    use yaml_rust::{YamlEmitter, YamlLoader};

    #[test]
    fn test_build_matrix() {
        let (u, nu, refs, reads, _) = build_matrix("tests/minimal_test.sam", None).unwrap();
        
        assert_eq!(refs.len(), 2, "Should have 2 references");
        assert_eq!(reads.len(), 3, "Should have 3 reads");
        assert_eq!(u.len(), 2, "Should have 2 unique reads (read1, read2)");
        assert_eq!(nu.len(), 1, "Should have 1 non-unique read (read3)");
        
        // Verify reference names
        assert!(refs.contains(&"ref1".to_string()));
        assert!(refs.contains(&"ref2".to_string()));
        
        // Verify read names
        assert!(reads.contains(&"read1".to_string()));
        assert!(reads.contains(&"read2".to_string()));
        assert!(reads.contains(&"read3".to_string()));
        
        // Verify data structure integrity
        for (_, (ref_idx, score)) in &u {
            assert!(*ref_idx >= 0 && *ref_idx < refs.len() as i32, "U matrix ref index should be valid");
            assert!(*score > 0.0, "U matrix scores should be positive");
        }
        
        for (_, (ref_indices, scores, normalized_scores, max_score)) in &nu {
            assert!(ref_indices.len() > 1, "NU entries should map to multiple references");
            assert_eq!(ref_indices.len(), scores.len(), "NU ref_indices and scores should have same length");
            assert_eq!(scores.len(), normalized_scores.len(), "NU scores and normalized_scores should have same length");
            
            // Verify normalized scores sum to approximately 1.0
            let sum: f64 = normalized_scores.iter().sum();
            assert!((sum - 1.0).abs() < 0.001, "Normalized scores should sum to 1.0");
            
            assert!(*max_score > 0.0, "Max score should be positive");
            
            // Verify all ref indices are valid
            for &ref_idx in ref_indices {
                assert!(ref_idx >= 0 && ref_idx < refs.len() as i32, "NU ref index should be valid");
            }
        }
        
        // Verify no duplicate references or reads
        let mut unique_refs = refs.clone();
        unique_refs.sort();
        unique_refs.dedup();
        assert_eq!(refs.len(), unique_refs.len(), "References should be unique");
        
        let mut unique_reads = reads.clone();
        unique_reads.sort();
        unique_reads.dedup();
        assert_eq!(reads.len(), unique_reads.len(), "Reads should be unique");
    }
}