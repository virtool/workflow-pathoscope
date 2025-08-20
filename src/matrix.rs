use crate::parse_sam::*;
use crate::{UniqueReads, MultiMappingReads};
use std::collections::HashMap;

pub fn build_matrix(
    alignment_path: &str,
    p_score_cutoff: Option<f64>,
) -> Result<
    (
        UniqueReads,
        MultiMappingReads,
        Vec<String>,
        Vec<String>,
        Vec<MinimalAlignment>,
    ),
    String,
> {
    let mut h_read_id: HashMap<String, i32> = HashMap::new();
    let mut h_ref_id: HashMap<String, i32> = HashMap::new();

    let mut refs: Vec<String> = Vec::new();
    let mut reads: Vec<String> = Vec::new();
    let mut minimal_alignments: Vec<MinimalAlignment> = Vec::new();

    let mut ref_count: i32 = 0;
    let mut read_count: i32 = 0;

    let mut max_score: f64 = 0.0;
    let mut min_score: f64 = 0.0;

    // First pass: collect all alignments per read
    let mut read_alignments: HashMap<i32, Vec<(i32, f64)>> = HashMap::new();

    // Parse all valid alignment lines from the file (SAM or BAM)
    let sam_lines = parse_alignment(alignment_path, p_score_cutoff)?;

    for new_line in sam_lines {

        // Store minimal alignment data for later use in coverage calculation
        let minimal_alignment = MinimalAlignment {
            read_idx: -1,  // Will be updated after read index is determined
            ref_idx: -1,   // Will be updated after ref index is determined
            position: new_line.position,
            read_length: new_line.read_length as u16,
        };
        minimal_alignments.push(minimal_alignment);

        let score = new_line.score.ok_or("Missing score in SAM line")?;
        min_score = score.min(min_score);
        max_score = score.max(max_score);

        // Get or create reference index
        let ref_index = *h_ref_id.get(&new_line.ref_id).unwrap_or(&-1);
        let ref_index = if ref_index == -1 {
            let new_ref_index = ref_count;
            h_ref_id.insert(new_line.ref_id.clone(), new_ref_index);
            refs.push(new_line.ref_id.clone());
            ref_count += 1;
            new_ref_index
        } else {
            ref_index
        };

        // Get or create read index
        let read_index = *h_read_id.get(&new_line.read_id).unwrap_or(&-1);
        let read_index = if read_index == -1 {
            let new_read_index = read_count;
            h_read_id.insert(new_line.read_id.clone(), new_read_index);
            reads.push(new_line.read_id.clone());
            read_count += 1;
            new_read_index
        } else {
            read_index
        };

        // Update the minimal alignment with the correct indices
        if let Some(last_alignment) = minimal_alignments.last_mut() {
            last_alignment.read_idx = read_index;
            last_alignment.ref_idx = ref_index;
        }

        // Add alignment to read's collection (skip duplicates)
        let alignments = read_alignments.entry(read_index).or_default();
        if !alignments.iter().any(|(ref_idx, _)| *ref_idx == ref_index) {
            alignments.push((ref_index, score));
        }
    }

    // Second pass: classify reads as unique or non-unique and build final data structures
    let mut u_temp: MultiMappingReads = std::collections::HashMap::new();
    let mut nu: MultiMappingReads = std::collections::HashMap::new();

    for (read_index, alignments) in read_alignments {
        if alignments.len() == 1 {
            // Unique read: maps to exactly one reference
            let (ref_index, score) = alignments[0];
            u_temp.insert(
                read_index,
                (vec![ref_index], vec![score], vec![score], score),
            );
        } else {
            // Non-unique read: maps to multiple references
            let ref_indices: Vec<i32> = alignments.iter().map(|(ref_idx, _)| *ref_idx).collect();
            let scores: Vec<f64> = alignments.iter().map(|(_, score)| *score).collect();
            let max_score = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            
            nu.insert(
                read_index,
                (ref_indices, scores, vec![], max_score),
            );
        }
    }

    let (u, mut nu) = rescale_samscore(u_temp, nu, max_score, min_score);

    let mut u_return: UniqueReads = std::collections::HashMap::new();

    for k in u.keys() {
        if let Some(u_entry) = u.get(k) {
            if let (Some(first_ref), Some(first_score)) = (u_entry.0.first(), u_entry.1.first()) {
                u_return.insert(*k, (*first_ref, *first_score));
            }
        }
    }

    let nu_keys: Vec<i32> = nu.keys().cloned().collect();
    for k in nu_keys {
        if let Some(nu_entry) = nu.get(&k) {
            let p_score_sum = nu_entry.1.iter().sum::<f64>();
            if let Some(nu_entry_mut) = nu.get_mut(&k) {
                nu_entry_mut.2 = nu_entry_mut
                    .1
                    .iter()
                    .map(|data| data / p_score_sum)
                    .collect();
            }
        }
    }

    Ok((u_return, nu, refs, reads, minimal_alignments))
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
        
        // With the score parsing fix, we should now process all 4 SAM lines correctly
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