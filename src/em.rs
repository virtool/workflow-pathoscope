use crate::{UniqueReads, MultiMappingReads};
use log::{info, warn};

pub fn compute_best_hit(
    u: &UniqueReads,
    nu: &MultiMappingReads,
    refs: &[String],
    reads: &[String],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ref_count = refs.len();
    let mut best_hit_reads = vec![0.0; ref_count];
    let mut level_1_reads = vec![0.0; ref_count];
    let mut level_2_reads = vec![0.0; ref_count];

    for i in u.keys() {
        if let Some(u_entry) = u.get(i) {
            let ref_idx = u_entry.0 as usize;
            if ref_idx < best_hit_reads.len() {
                best_hit_reads[ref_idx] += 1.0;
            }
            if ref_idx < level_1_reads.len() {
                level_1_reads[ref_idx] += 1.0;
            }
        }
    }

    for i in nu.keys() {
        if let Some(z) = nu.get(i) {
            let ind = &z.0;
            let x_norm = &z.2;
            let best_ref = x_norm.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let mut num_best_ref = 0;

            for score in x_norm.iter() {
                if *score == best_ref {
                    num_best_ref += 1;
                }
            }

            num_best_ref = match num_best_ref {
                0 => 1,
                _ => num_best_ref,
            };

            for (j, score) in x_norm.iter().enumerate() {
                if let Some(&ref_idx) = ind.get(j) {
                    let ref_idx = ref_idx as usize;
                    
                    if *score == best_ref
                        && ref_idx < best_hit_reads.len() {
                            best_hit_reads[ref_idx] += 1.0 / num_best_ref as f64;
                        }

                    if *score >= 0.5 && ref_idx < level_1_reads.len() {
                        level_1_reads[ref_idx] += 1.0;
                    } else if *score >= 0.01 && ref_idx < level_2_reads.len() {
                        level_2_reads[ref_idx] += 1.0;
                    }
                }
            }
        }
    }

    let read_count = reads.len();

    let best_hit: Vec<f64> = best_hit_reads
        .iter()
        .map(|val| *val / read_count as f64)
        .collect();
    let level1: Vec<f64> = level_1_reads
        .iter()
        .map(|val| *val / read_count as f64)
        .collect();
    let level2: Vec<f64> = level_2_reads
        .iter()
        .map(|val| *val / read_count as f64)
        .collect();

    (best_hit_reads, best_hit, level1, level2)
}

/// Parameters used throughout the EM algorithm iterations
struct EMParameters {
    pi: Vec<f64>,
    theta: Vec<f64>,
    pi_sum_0: Vec<f64>,
    u_total: f64,
    nu_total: f64,
    prior_weight: f64,
    nu_length: usize,
}

/// Initialize parameters for the EM algorithm
/// 
/// Computes initial pi and theta values (uniform distribution), calculates weight
/// statistics from unique and multi-mapping reads, and sets up data structures
/// needed for EM iterations.
/// 
/// # Arguments
/// * `u` - Unique read mappings (read_id -> (ref_id, score))
/// * `nu` - Multi-mapping read data (read_id -> (ref_ids, scores, normalized_scores, max_score))
/// * `genome_count` - Number of reference genomes
/// 
/// # Returns
/// EMParameters struct containing all initialization data
fn initialize_em_parameters(
    u: &UniqueReads,
    nu: &MultiMappingReads,
    genome_count: usize,
) -> EMParameters {
    let pi = vec![1.0 / genome_count as f64; genome_count];
    let theta = pi.clone();
    let mut pi_sum_0 = vec![0.0; genome_count];

    let u_weights: Vec<f64> = u.iter().map(|entry| (entry.1).1).collect();
    let mut max_u_weights = 0.0;
    let mut u_total = 0.0;

    if !u_weights.is_empty() {
        max_u_weights = u_weights.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        u_total = u_weights.iter().sum();
    }

    for i in u.keys() {
        if let Some(u_entry) = u.get(i) {
            let ref_idx = u_entry.0 as usize;
            if ref_idx < pi_sum_0.len() {
                pi_sum_0[ref_idx] += u_entry.1;
            }
        }
    }

    let nu_weights: Vec<f64> = nu.iter().map(|entry| (entry.1).3).collect();
    let mut max_nu_weights = 0.0;
    let mut nu_total = 0.0;

    if !nu_weights.is_empty() {
        max_nu_weights = nu_weights.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        nu_total = nu_weights.iter().sum();
    }

    let prior_weight = f64::max(max_u_weights, max_nu_weights);
    let nu_length = if nu.is_empty() { 1 } else { nu.len() };

    EMParameters {
        pi,
        theta,
        pi_sum_0,
        u_total,
        nu_total,
        prior_weight,
        nu_length,
    }
}

/// Expectation step of the EM algorithm
/// 
/// Updates posterior probabilities for multi-mapping reads based on current
/// pi and theta estimates. For each multi-mapping read, calculates the 
/// probability of assignment to each mapped reference.
/// 
/// # Arguments
/// * `nu` - Multi-mapping read data (modified in place with updated probabilities)
/// * `pi` - Current genome abundance estimates
/// * `theta` - Current multi-mapping probability estimates
/// * `theta_sum` - Output buffer for accumulating weighted assignments (modified in place)
fn em_e_step(
    nu: &mut MultiMappingReads,
    pi: &[f64],
    theta: &[f64],
    theta_sum: &mut [f64],
) {
    let nu_keys: Vec<i32> = nu.keys().cloned().collect();
    for j in nu_keys {
        if let Some(z) = nu.get(&j).cloned() {
            //A set of any genome mapping with j
            let ind = &z.0;

            //Get relevant pis for the read
            let pi_temp: Vec<f64> = ind
                .iter()
                .filter_map(|&val| pi.get(val as usize).cloned())
                .collect();

            //Get relevant thetas for the read
            let theta_temp: Vec<f64> = ind
                .iter()
                .filter_map(|&val| theta.get(val as usize).cloned())
                .collect();

            //Calculate non-normalized xs
            let mut x_temp: Vec<f64> = Vec::with_capacity(ind.len());

            for k in 0..ind.len().min(pi_temp.len()).min(theta_temp.len()) {
                if let Some(&score) = z.1.get(k) {
                    x_temp.push(pi_temp[k] * theta_temp[k] * score);
                }
            }

            let x_sum: f64 = x_temp.iter().sum();

            //Avoid dividing by 0 at all times
            let x_norm: Vec<f64> = if x_sum == 0.0 {
                vec![0.0; x_temp.len()]
            } else {
                x_temp.iter().map(|val| val / x_sum).collect()
            };

            //Update x in nu
            if let Some(nu_entry) = nu.get_mut(&j) {
                nu_entry.2.clone_from(&x_norm);

                // Only update theta_sum if we have meaningful scores (optimization)
                if x_sum > 0.0 {
                    for (k, &ref_idx) in ind.iter().enumerate() {
                        if let Some(&x_val) = x_norm.get(k) {
                            let idx = ref_idx as usize;
                            if idx < theta_sum.len() {
                                theta_sum[idx] += x_val * nu_entry.3;
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Maximization step of the EM algorithm
/// 
/// Updates pi and theta parameters to maximize the expected log-likelihood
/// based on the current assignment probabilities from the E-step.
/// 
/// # Arguments
/// * `params` - EM parameters containing priors and weights
/// * `theta_sum` - Accumulated weighted assignments from E-step
/// * `pi_prior` - Prior weight for pi estimation (Dirichlet prior)
/// * `theta_prior` - Prior weight for theta estimation
/// 
/// # Returns
/// Tuple of (new_pi, new_theta) vectors with updated parameter estimates
fn em_m_step(
    params: &EMParameters,
    theta_sum: &[f64],
    pi_prior: f64,
    theta_prior: f64,
) -> (Vec<f64>, Vec<f64>) {
    let pi_sum: Vec<f64> = theta_sum
        .iter()
        .enumerate()
        .map(|(idx, _)| theta_sum[idx] + params.pi_sum_0[idx])
        .collect();
    let pip = pi_prior * params.prior_weight;

    //update pi
    let new_pi = pi_sum
        .iter()
        .map(|val| ((*val) + pip) / (params.u_total + params.nu_total + (pip * pi_sum.len() as f64)))
        .collect();

    let theta_p = theta_prior * params.prior_weight;

    let nu_total_div = if params.nu_total == 0.0 { 1.0 } else { params.nu_total };

    let new_theta = theta_sum
        .iter()
        .map(|val| (*val + theta_p) / (nu_total_div + (theta_p * theta_sum.len() as f64)))
        .collect();

    (new_pi, new_theta)
}

/// Calculate convergence metric as L1 norm between old and new pi vectors
/// 
/// # Arguments
/// * `pi_old` - Previous iteration's pi values
/// * `pi_new` - Current iteration's pi values
/// 
/// # Returns
/// Sum of absolute differences between old and new pi values
fn calculate_convergence(pi_old: &[f64], pi_new: &[f64]) -> f64 {
    pi_old.iter()
        .zip(pi_new.iter())
        .map(|(old, new)| (old - new).abs())
        .sum()
}

/// Check if EM algorithm has converged and log progress
/// 
/// Determines convergence based on the change in pi values between iterations.
/// Also handles progress logging and divergence detection.
/// 
/// # Arguments
/// * `iteration` - Current iteration number (0-indexed)
/// * `cutoff` - Current convergence metric value
/// * `epsilon` - Convergence threshold
/// * `nu_length` - Number of multi-mapping reads
/// 
/// # Returns
/// true if converged, false otherwise
fn check_convergence(
    iteration: i32,
    cutoff: f64,
    epsilon: f64,
    nu_length: usize,
) -> bool {
    // Log convergence progress
    if iteration == 0 || iteration % 10 == 9 || cutoff <= epsilon {
        info!("em iteration {}: convergence delta = {:.2e}", iteration + 1, cutoff);
    }

    if cutoff <= epsilon || nu_length == 1 {
        info!("em converged after {} iterations (delta: {:.2e})", iteration + 1, cutoff);
        return true;
    }

    // Detect potential divergence
    if iteration > 10 && cutoff > 1e-2 {
        info!("em may be diverging at iteration {} (delta: {:.2e})", iteration + 1, cutoff);
    }

    false
}

pub fn em(
    u: &UniqueReads,
    mut nu: MultiMappingReads,
    genomes: &[String],
    max_iter: i32,
    epsilon: f64,
    pi_prior: f64,
    theta_prior: f64,
) -> (
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    MultiMappingReads,
) {
    let genome_count = genomes.len();
    let params = initialize_em_parameters(u, &nu, genome_count);
    let mut pi = params.pi.clone();
    let mut theta = params.theta.clone();
    let mut init_pi = Vec::new();

    //EM iterations
    for i in 0..max_iter {
        let pi_old = pi.clone();
        let mut theta_sum = vec![0.0; genome_count];

        //E step
        em_e_step(&mut nu, &pi, &theta, &mut theta_sum);

        //M step
        let (new_pi, new_theta) = em_m_step(
            &params,
            &theta_sum,
            pi_prior,
            theta_prior,
        );
        
        pi = new_pi;
        theta = new_theta;

        if i == 0 {
            init_pi.clone_from(&pi);
        }

        let cutoff = calculate_convergence(&pi_old, &pi);
        if check_convergence(i, cutoff, epsilon, params.nu_length) {
            break;
        }
    }

    (init_pi, pi, theta, nu)
}

pub fn find_updated_score(
    nu: &MultiMappingReads,
    read_index: i32,
    ref_index: i32,
) -> f64 {
    let v1 = match nu.get(&read_index) {
        Some(val) => val,
        None => return 0.0,
    };

    // Find the index of ref_index in the reference list
    for (i, &el) in v1.0.iter().enumerate() {
        if el == ref_index {
            return v1.2.get(i).copied().unwrap_or(0.0);
        }
    }
    
    // If ref_index not found, return 0.0
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_updated_score() {
        let mut nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Read 0 maps to refs [0, 1] with normalized scores [0.7, 0.3]
        nu.insert(0, (
            vec![0, 1], 
            vec![80.0, 20.0], 
            vec![0.7, 0.3], 
            100.0
        ));
        
        // Read 1 maps to refs [1, 2] with normalized scores [0.6, 0.4]
        nu.insert(1, (
            vec![1, 2], 
            vec![60.0, 40.0], 
            vec![0.6, 0.4], 
            100.0
        ));

        // Test valid lookups
        assert_eq!(find_updated_score(&nu, 0, 0), 0.7, "Read 0, ref 0 should return 0.7");
        assert_eq!(find_updated_score(&nu, 0, 1), 0.3, "Read 0, ref 1 should return 0.3");
        assert_eq!(find_updated_score(&nu, 1, 1), 0.6, "Read 1, ref 1 should return 0.6");
        assert_eq!(find_updated_score(&nu, 1, 2), 0.4, "Read 1, ref 2 should return 0.4");

        // Test invalid read index
        assert_eq!(find_updated_score(&nu, 999, 0), 0.0, "Non-existent read should return 0.0");

        // Test invalid ref index for existing read
        assert_eq!(find_updated_score(&nu, 0, 999), 0.0, "Non-existent ref should return 0.0");
        assert_eq!(find_updated_score(&nu, 0, 2), 0.0, "Read 0 doesn't map to ref 2, should return 0.0");
    }

    #[test]
    fn test_find_updated_score_empty_data() {
        let nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Test with empty HashMap
        assert_eq!(find_updated_score(&nu, 0, 0), 0.0, "Empty nu should return 0.0");
    }

    #[test]
    fn test_find_updated_score_edge_cases() {
        let mut nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Read with empty ref vectors
        nu.insert(0, (vec![], vec![], vec![], 0.0));
        assert_eq!(find_updated_score(&nu, 0, 0), 0.0, "Read with empty refs should return 0.0");
        
        // Read with mismatched vector lengths (should not happen in practice, but test robustness)
        nu.insert(1, (
            vec![0], 
            vec![100.0], 
            vec![], // Empty normalized scores
            100.0
        ));
        assert_eq!(find_updated_score(&nu, 1, 0), 0.0, "Empty normalized scores should return 0.0");
    }

    #[test]
    fn test_find_updated_score_ref_not_found() {
        let mut nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Read 0 maps to refs [1, 2] with scores [0.8, 0.2]
        nu.insert(0, (
            vec![1, 2], 
            vec![80.0, 20.0], 
            vec![0.8, 0.2], 
            100.0
        ));

        // Looking for ref 0 (which doesn't exist) should return 0.0
        assert_eq!(find_updated_score(&nu, 0, 0), 0.0, "Should return 0.0 when ref not found");
    }

    #[test]
    fn test_em_basic_convergence() {
        // Simple test case: 2 references, 2 reads
        let mut u: UniqueReads = std::collections::HashMap::new();
        let mut nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Read 0 uniquely maps to ref 0 with score 100.0
        u.insert(0, (0, 100.0));
        
        // Read 1 multi-maps to both refs with equal scores
        nu.insert(1, (
            vec![0, 1],           // Maps to refs 0 and 1
            vec![50.0, 50.0],     // Equal raw scores
            vec![0.5, 0.5],       // Equal normalized scores (will be updated by EM)
            100.0                 // Total score
        ));
        
        let genomes = vec!["ref0".to_string(), "ref1".to_string()];
        
        // Run EM algorithm
        let (init_pi, final_pi, theta, updated_nu) = em(
            &u, 
            nu, 
            &genomes, 
            10,    // max_iter
            1e-6,  // epsilon
            0.0,   // pi_prior
            0.0    // theta_prior
        );
        
        // Test return value structure
        assert_eq!(init_pi.len(), 2, "init_pi should have 2 elements");
        assert_eq!(final_pi.len(), 2, "final_pi should have 2 elements");
        assert_eq!(theta.len(), 2, "theta should have 2 elements");
        
        // Test that pi values sum to approximately 1.0
        let pi_sum: f64 = final_pi.iter().sum();
        assert!((pi_sum - 1.0).abs() < 0.01, "Pi values should sum to ~1.0, got {}", pi_sum);
        
        // Test that theta values sum to approximately 1.0
        let theta_sum: f64 = theta.iter().sum();
        assert!((theta_sum - 1.0).abs() < 0.01, "Theta values should sum to ~1.0, got {}", theta_sum);
        
        // Since read 0 uniquely maps to ref 0, ref 0 should have higher pi
        assert!(final_pi[0] > final_pi[1], "ref0 should have higher pi due to unique read, pi0={}, pi1={}", final_pi[0], final_pi[1]);
        
        // All values should be positive
        assert!(final_pi[0] > 0.0 && final_pi[1] > 0.0, "All pi values should be positive");
        assert!(theta[0] > 0.0 && theta[1] > 0.0, "All theta values should be positive");
        
        // Verify the multi-mapping read's scores were updated
        let updated_read1 = updated_nu.get(&1).expect("Read 1 should exist in updated nu");
        let updated_scores = &updated_read1.2;
        assert_eq!(updated_scores.len(), 2, "Updated scores should have 2 elements");
        
        // The sum of normalized scores should be 1.0
        let score_sum: f64 = updated_scores.iter().sum();
        assert!((score_sum - 1.0).abs() < 1e-10, "Normalized scores should sum to 1.0, got {}", score_sum);
    }

    #[test]
    fn test_em_empty_input() {
        // Test with completely empty input
        let u: UniqueReads = std::collections::HashMap::new();
        let nu: MultiMappingReads = std::collections::HashMap::new();
        let genomes = vec!["ref0".to_string(), "ref1".to_string()];
        
        let (init_pi, final_pi, theta, updated_nu) = em(
            &u, 
            nu, 
            &genomes, 
            10, 1e-6, 0.0, 0.0
        );
        
        // Should return equal probabilities for all references
        assert_eq!(init_pi.len(), 2);
        assert_eq!(final_pi.len(), 2);
        assert_eq!(theta.len(), 2);
        
        // With no data, the algorithm produces NaN values due to division by zero
        // This documents the current behavior - empty input results in NaN
        assert!(init_pi[0].is_nan(), "Empty input produces NaN for init_pi[0]");
        assert!(init_pi[1].is_nan(), "Empty input produces NaN for init_pi[1]");
        assert!(final_pi[0].is_nan(), "Empty input produces NaN for final_pi[0]");
        assert!(final_pi[1].is_nan(), "Empty input produces NaN for final_pi[1]");
        
        // Nu should remain empty
        assert!(updated_nu.is_empty(), "Updated nu should remain empty");
    }

    #[test]
    fn test_em_only_unique_reads() {
        // Test with only unique reads (no multi-mapping)
        let mut u: UniqueReads = std::collections::HashMap::new();
        let nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Add unique reads: 2 to ref0, 1 to ref1
        u.insert(0, (0, 100.0));  // read 0 -> ref 0
        u.insert(1, (0, 100.0));  // read 1 -> ref 0  
        u.insert(2, (1, 100.0));  // read 2 -> ref 1
        
        let genomes = vec!["ref0".to_string(), "ref1".to_string()];
        
        let (init_pi, final_pi, theta, updated_nu) = em(
            &u, 
            nu, 
            &genomes, 
            10, 1e-6, 0.0, 0.0
        );
        
        // ref0 should have higher pi due to more reads (2 vs 1)
        assert!(final_pi[0] > final_pi[1], "ref0 should have higher pi with more unique reads");
        
        // Pi should reflect the read distribution: ref0 gets 2/3, ref1 gets 1/3
        assert!((final_pi[0] - 2.0/3.0).abs() < 0.01, "ref0 should get ~2/3 of probability");
        assert!((final_pi[1] - 1.0/3.0).abs() < 0.01, "ref1 should get ~1/3 of probability");
        
        // Nu should remain empty
        assert!(updated_nu.is_empty(), "Nu should remain empty with only unique reads");
    }

    #[test]
    fn test_em_only_multi_mapping_reads() {
        // Test with only multi-mapping reads (no unique reads)
        let u: UniqueReads = std::collections::HashMap::new();
        let mut nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Add multi-mapping reads
        nu.insert(0, (vec![0, 1], vec![60.0, 40.0], vec![0.6, 0.4], 100.0));
        nu.insert(1, (vec![0, 1], vec![30.0, 70.0], vec![0.3, 0.7], 100.0));
        
        let genomes = vec!["ref0".to_string(), "ref1".to_string()];
        
        let (init_pi, final_pi, theta, updated_nu) = em(
            &u, 
            nu, 
            &genomes, 
            10, 1e-6, 0.0, 0.0
        );
        
        // Should converge to some distribution
        let pi_sum: f64 = final_pi.iter().sum();
        assert!((pi_sum - 1.0).abs() < 0.01, "Pi should sum to 1.0");
        
        // All probabilities should be positive
        assert!(final_pi[0] > 0.0 && final_pi[1] > 0.0, "All pi values should be positive");
        assert!(theta[0] > 0.0 && theta[1] > 0.0, "All theta values should be positive");
        
        // Updated nu should have normalized scores
        for (_, entry) in updated_nu.iter() {
            let score_sum: f64 = entry.2.iter().sum();
            assert!((score_sum - 1.0).abs() < 1e-10, "Each read's scores should sum to 1.0");
        }
    }

    #[test]
    fn test_em_integration_with_real_sam_data() {
        use crate::build_matrix;
        
        // Use real SAM data with multi-mapping reads to test the full pipeline
        let sam_path = "example/rust/test_em_with_multimapping.sam";
        
        // Build matrix from SAM file
        let (u, nu, refs, reads, _minimal_alignments) = build_matrix(sam_path, None)
            .expect("Failed to build matrix from test SAM file");
        
        // Run EM algorithm
        let (init_pi, final_pi, theta, updated_nu) = em(
            &u, 
            nu, 
            &refs, 
            50,    // max_iter
            1e-7,  // epsilon  
            0.0,   // pi_prior
            0.0    // theta_prior
        );
        
        // Test return value structure
        assert_eq!(init_pi.len(), refs.len(), "init_pi should match number of refs");
        assert_eq!(final_pi.len(), refs.len(), "final_pi should match number of refs");
        assert_eq!(theta.len(), refs.len(), "theta should match number of refs");
        
        // Test normalization
        let pi_sum: f64 = final_pi.iter().sum();
        assert!((pi_sum - 1.0).abs() < 0.01, "Pi values should sum to ~1.0, got {}", pi_sum);
        
        let theta_sum: f64 = theta.iter().sum();
        // Theta represents multi-mapping read probabilities, should sum to ~1.0 when we have multi-mapping reads
        if updated_nu.len() > 0 {
            assert!((theta_sum - 1.0).abs() < 0.01, "Theta values should sum to ~1.0 when multi-mapping reads exist, got {}", theta_sum);
        } else {
            // If no multi-mapping reads, theta can be 0 or very small
            assert!(theta_sum >= 0.0, "Theta sum should be non-negative, got {}", theta_sum);
        }
        
        // All probabilities should be positive
        for (i, &pi_val) in final_pi.iter().enumerate() {
            assert!(pi_val > 0.0, "Pi[{}] should be positive, got {}", i, pi_val);
        }
        
        for (i, &theta_val) in theta.iter().enumerate() {
            assert!(theta_val >= 0.0, "Theta[{}] should be non-negative, got {}", i, theta_val);
        }
        
        // Test that multi-mapping reads have normalized scores
        for (read_id, entry) in updated_nu.iter() {
            let score_sum: f64 = entry.2.iter().sum();
            assert!((score_sum - 1.0).abs() < 1e-10, 
                "Read {} scores should sum to 1.0, got {}", read_id, score_sum);
            
            // All individual scores should be non-negative
            for (j, &score) in entry.2.iter().enumerate() {
                assert!(score >= 0.0, 
                    "Read {} score[{}] should be non-negative, got {}", read_id, j, score);
            }
        }
        
        // Test that the number of references and reads matches expectations
        assert!(refs.len() >= 1, "Should have at least one reference");
        assert!(reads.len() >= 1, "Should have at least one read");
        
        // Test that we have both unique and multi-mapping reads
        assert!(u.len() > 0, "Should have some unique reads, got {}", u.len());
        assert!(updated_nu.len() > 0, "Should have some multi-mapping reads, got {}", updated_nu.len());
        
        // Test that theta is meaningful (positive) since we have multi-mapping reads
        let theta_sum: f64 = theta.iter().sum();
        assert!(theta_sum > 0.0, "Theta should be positive when multi-mapping reads exist, got {}", theta_sum);
        
        // Test specific expectations for this SAM file:
        // - Should have 2 references (ref1, ref2) 
        // - Should have 6 total reads (4 unique + 2 multi-mapping)
        assert_eq!(refs.len(), 2, "Should have 2 references");
        assert_eq!(reads.len(), 6, "Should have 6 total reads");
        
        // Verify that ref1 should have highest pi (gets 3 unique reads vs 1 for ref2)
        let max_pi_idx = final_pi.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
        assert_eq!(max_pi_idx, 0, "ref1 should have highest pi due to more unique reads");
        
        println!("Integration test completed successfully:");
        println!("  Refs: {}, Reads: {}", refs.len(), reads.len());
        println!("  Unique reads: {}, Multi-mapping reads: {}", u.len(), updated_nu.len());
        println!("  Final pi: {:?}", final_pi);
        println!("  Final theta: {:?}", theta);
    }

    #[test]
    fn test_em_with_zero_score_reads() {
        // Test with multi-mapping reads that have zero scores
        let u: UniqueReads = std::collections::HashMap::new();
        let mut nu: MultiMappingReads = std::collections::HashMap::new();
        
        // Add multi-mapping reads with zero scores
        nu.insert(0, (vec![0, 1], vec![0.0, 0.0], vec![0.0, 0.0], 0.0)); // All zero scores
        nu.insert(1, (vec![0, 1], vec![50.0, 30.0], vec![0.625, 0.375], 80.0)); // Normal scores
        nu.insert(2, (vec![0, 1], vec![0.0, 100.0], vec![0.0, 1.0], 100.0)); // One zero score
        
        let genomes = vec!["ref0".to_string(), "ref1".to_string()];
        
        let (init_pi, final_pi, theta, updated_nu) = em(
            &u, 
            nu, 
            &genomes, 
            10, 1e-6, 0.0, 0.0
        );
        
        // Should still converge despite zero scores
        let pi_sum: f64 = final_pi.iter().sum();
        assert!((pi_sum - 1.0).abs() < 0.01, "Pi should sum to 1.0 even with zero scores");
        
        // Check that reads with all zero scores don't break the algorithm
        assert!(final_pi[0] >= 0.0 && final_pi[1] >= 0.0, "All pi values should be non-negative");
        assert!(theta[0] >= 0.0 && theta[1] >= 0.0, "All theta values should be non-negative");
        
        // Verify that updated_nu contains the reads
        assert_eq!(updated_nu.len(), 3, "Should have 3 multi-mapping reads");
        
        // Check that reads with meaningful scores still have normalized scores
        for (read_id, entry) in updated_nu.iter() {
            if *read_id == 1 || *read_id == 2 { // Reads with non-zero scores
                let score_sum: f64 = entry.2.iter().sum();
                assert!((score_sum - 1.0).abs() < 1e-10, 
                    "Read {} with non-zero scores should sum to 1.0, got {}", read_id, score_sum);
            }
        }
        
        // Read 0 with all zero scores should have specific behavior
        if let Some(zero_read) = updated_nu.get(&0) {
            let score_sum: f64 = zero_read.2.iter().sum();
            // Current behavior: zero scores get normalized to [0.0, 0.0] which sums to 0.0
            assert_eq!(score_sum, 0.0, "Read with all zero scores should sum to 0.0");
        }
    }
}
