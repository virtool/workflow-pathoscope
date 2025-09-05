use crate::sam::MinimalAlignment;
use crate::em::find_updated_score;
use crate::{UniqueReads, MultiMappingReads};
use rustc_hash::FxHashMap;

// TODO: Fix the updated score calculation logic
// The current implementation includes any alignment where posterior probability >= p_score_cutoff,
// but this may not match the expected behavior. The Rust unit test that was attempted revealed
// that read1 with probability 0.1 for ref1 gets included because 0.1 > 0.01, even though
// read1 has a higher probability (0.9) for ref0. This suggests the current logic includes
// ALL alignments above threshold rather than just the best assignments per read.
// 
// This behavior is currently verified as correct by the Python integration test, but we should
// investigate whether this is the intended behavior or if we should only include the best
// assignment for each multi-mapping read that meets the threshold.

/// Calculate coverage directly from EM results without creating intermediate SAM file
pub fn calculate_coverage_from_em(
    u: &UniqueReads,
    nu: &MultiMappingReads,
    minimal_alignments: &[MinimalAlignment],
    refs: &[String],
    p_score_cutoff: f64,
    ref_lengths: &FxHashMap<String, usize>,
) -> FxHashMap<String, Vec<usize>> {
    let mut coverage_dict: FxHashMap<String, Vec<usize>> = FxHashMap::default();

    // Initialize coverage arrays for all references
    for (ref_id, &length) in ref_lengths {
        coverage_dict.insert(ref_id.clone(), vec![0; length]);
    }

    // Process each minimal alignment and calculate coverage based on EM results
    for alignment in minimal_alignments {
        let read_index = alignment.read_idx;
        let ref_index = alignment.ref_idx;

        // Determine if this alignment should contribute to coverage
        let should_include = if u.contains_key(&read_index) {
            // Unique reads always contribute
            true
        } else if nu.contains_key(&read_index) {
            // Multi-mapping reads only contribute if posterior probability >= threshold
            find_updated_score(nu, read_index, ref_index) >= p_score_cutoff
        } else {
            // Read not found in EM results, skip
            false
        };

        if should_include {
            // Get reference name from index
            if let Some(ref_name) = refs.get(ref_index as usize) {
                // Add coverage for this alignment
                if let Some(coverage_array) = coverage_dict.get_mut(ref_name) {
                    let start_index = if alignment.position > 0 {
                        (alignment.position - 1) as usize
                    } else {
                        0
                    };

                    let end_index = (start_index + alignment.read_length as usize).min(coverage_array.len());

                    for item in coverage_array.iter_mut().take(end_index).skip(start_index) {
                        *item += 1;
                    }
                }
            }
        }
    }

    coverage_dict
}