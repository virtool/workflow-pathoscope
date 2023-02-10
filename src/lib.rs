use eliminate_subtraction::{
    check_should_eliminate, find_sam_align_score, parse_subtraction_sam, read_lines,
};
use pyo3::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::Write,
};

#[pymodule]
///pyo3 interface
fn rust_utils(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_expectation_maximization, m)?)?;
    m.add_function(wrap_pyfunction!(run_eliminate_subtraction, m)?)?;
    return Ok(());
}

#[pyfunction]
///Entry point for eliminate_subtraction
pub fn run_eliminate_subtraction(
    _py: Python,
    isolate_sam_path: String,
    subtraction_sam_path: String,
    output_sam_path: String,
) {
    let subtraction_scores = parse_subtraction_sam(&subtraction_sam_path);

    if let Ok(lines) = read_lines(isolate_sam_path) {
        let mut sam_file = File::create(output_sam_path).unwrap();
        let mut subtracted_read_ids: HashSet<String> = HashSet::new();

        for line in lines {
            if let Ok(l) = line {
                match l.chars().next() {
                    Some(c) => {
                        if c == '@' || c == '#' {
                            continue;
                        }
                    }
                    None => continue,
                };

                let first = l.chars().next().unwrap();

                if first == '@' || first == '#' {
                    continue;
                }

                let fields: Vec<&str> = l.split("\t").collect();

                if fields[2] == "*" {
                    continue;
                }

                let score = find_sam_align_score(&fields);

                let eliminate = check_should_eliminate(&subtraction_scores, &fields[0], score);

                if eliminate {
                    subtracted_read_ids.insert(fields[0].to_string());
                } else {
                    writeln!(&mut sam_file, "{}", l).unwrap();
                }
            }
        }

        let mut subtracted_read_ids_file = File::create("subtracted_read_ids.txt").unwrap();

        for read_id in subtracted_read_ids {
            writeln!(&mut subtracted_read_ids_file, "{}", read_id).unwrap();
        }
    }
}

mod eliminate_subtraction {
    use std::{
        collections::HashMap,
        fs::File,
        io::{self, BufRead},
        path::Path,
    };

    // Check if the passed read_id should be eliminated if its isolate score is
    // higher than the subtraction score.
    pub fn check_should_eliminate(
        subtraction_scores: &HashMap<String, f32>,
        read_id: &str,
        score: f32,
    ) -> bool {
        match subtraction_scores.get(read_id) {
            Some(subtraction_score) => &subtraction_score >= &&score,
            None => false,
        }
    }

    /// Find the Pathoscope alignment score for a SAM line.
    ///
    /// # Arguments
    /// * `fields` - The SAM fields as a vector.
    ///
    pub fn find_sam_align_score(fields: &Vec<&str>) -> f32 {
        let read_length = fields[9].chars().count() as f32;
        let mut a_score: f32 = 0.0;

        for field in fields {
            if field.starts_with("AS:i:") {
                a_score = field[5..].parse().unwrap();
                break;
            }
        }

        return a_score + read_length;
    }

    pub fn parse_subtraction_sam(path: &str) -> HashMap<String, f32> {
        let mut high_scores: HashMap<String, f32> = HashMap::new();

        if let Ok(lines) = read_lines(path) {
            for line in lines {
                if let Ok(l) = line {
                    let first = l.chars().next().unwrap();

                    if first == '@' || first == '#' {
                        continue;
                    }

                    let fields: Vec<&str> = l.split("\t").collect();

                    if fields[2] == "*" {
                        continue;
                    }

                    let score = find_sam_align_score(&fields);
                    high_scores.insert(fields[0].to_string(), score);
                }
            }
        }

        return high_scores;
    }

    pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        Ok(io::BufReader::new(file).lines())
    }
}

#[pyfunction]
///Entry point for expectation_maximization
pub fn run_expectation_maximization(
    _py: Python,
    sam_path: String,
    reassigned_path: String,
    p_score_cutoff: f64,
) -> (
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<String>,
    Vec<String>,
) {
    let (u, nu, refs, reads) = build_matrix::build_matrix(sam_path.as_str(), None);

    let (best_hit_initial_reads, best_hit_initial, level_1_initial, level_2_initial) =
        compute_best_hit(&u, &nu, &refs, &reads);

    let (init_pi, pi, _, nu) = em(&u, nu, &refs, 50, 1e-7, 0.0, 0.0);

    let (best_hit_final_reads, best_hit_final, level_1_final, level_2_final) =
        compute_best_hit(&u, &nu, &refs, &reads);

    rewrite_align::rewrite_align(
        &u,
        &nu,
        sam_path.as_str(),
        &p_score_cutoff,
        &reassigned_path,
    );

    return (
        best_hit_initial_reads,
        best_hit_initial,
        level_1_initial,
        level_2_initial,
        best_hit_final_reads,
        best_hit_final,
        level_1_final,
        level_2_final,
        init_pi,
        pi,
        refs,
        reads,
    );
}

mod build_matrix {
    use super::parse_sam::*;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;

    pub fn build_matrix(
        sam_path: &str,
        p_score_cutoff: Option<f64>,
    ) -> (
        HashMap<i32, (i32, f64)>,
        HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        Vec<String>,
        Vec<String>,
    ) {
        let mut u: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)> = HashMap::new();
        let mut nu: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)> = HashMap::new();

        let mut h_read_id: HashMap<String, i32> = HashMap::new();
        let mut h_ref_id: HashMap<String, i32> = HashMap::new();

        let mut refs: Vec<String> = Vec::new();
        let mut reads: Vec<String> = Vec::new();

        let mut ref_count: i32 = 0;
        let mut read_count: i32 = 0;

        let mut max_score: f64 = 0.0;
        let mut min_score: f64 = 0.0;

        let sam_file = File::open(sam_path).expect("Invalid file");
        let mut sam_reader = BufReader::new(sam_file);

        loop {
            let new_line: SamLine;

            match parse_sam(&mut sam_reader, p_score_cutoff) {
                ParseResult::Ok(data) => {
                    if data.score < p_score_cutoff {
                        continue;
                    } else {
                        new_line = data;
                    }
                }
                ParseResult::EOF => break,
                ParseResult::Err(msg) => {
                    println!("{}", msg);
                    panic!();
                }
                ParseResult::Ignore => continue,
            }

            min_score = new_line.score.unwrap().min(min_score);
            max_score = new_line.score.unwrap().max(max_score);

            let mut ref_index = *(h_ref_id.get(&new_line.ref_id).unwrap_or(&-1));

            if ref_index == -1 {
                ref_index = ref_count;
                h_ref_id.insert(new_line.ref_id.clone(), ref_index);
                refs.push(new_line.ref_id.clone());
                ref_count += 1;
            }

            let mut read_index = *(h_read_id.get(&new_line.read_id).unwrap_or(&-1));

            if read_index == -1 {
                read_index = read_count;
                h_read_id.insert(new_line.read_id.clone(), read_index);
                reads.push(new_line.read_id.clone());
                read_count += 1;

                u.insert(
                    read_index,
                    (
                        vec![ref_index],
                        vec![new_line.score.clone().unwrap()],
                        vec![new_line.score.clone().unwrap() as f64],
                        new_line.score.clone().unwrap(),
                    ),
                );
            } else {
                if u.contains_key(&read_index) {
                    if u.get(&read_index).unwrap().0.contains(&ref_index) {
                        continue;
                    }
                    nu.insert(read_index, u.get(&read_index).unwrap().clone());
                    u.remove(&read_index);
                }

                if nu.get(&read_index).unwrap().0.contains(&ref_index) {
                    continue;
                }

                nu.get_mut(&read_index).unwrap().0.push(ref_index);
                nu.get_mut(&read_index)
                    .unwrap()
                    .1
                    .push(new_line.score.unwrap());

                if new_line.score.unwrap() > nu.get(&read_index).unwrap().3 {
                    nu.get_mut(&read_index).unwrap().3 = new_line.score.unwrap();
                }
            }
        }

        let (u, mut nu) = rescale_samscore(u, nu, max_score, min_score);

        let mut u_return: HashMap<i32, (i32, f64)> = HashMap::new();

        for k in u.keys() {
            u_return.insert(
                *k,
                (
                    u.get(k).unwrap().0.get(0).unwrap().clone(),
                    u.get(k).unwrap().1.get(0).unwrap().clone(),
                ),
            );
        }

        for k in nu.clone().keys() {
            let p_score_sum = nu.get(k).unwrap().1.iter().sum::<f64>();

            nu.get_mut(k).unwrap().2 = nu
                .get(k)
                .unwrap()
                .1
                .iter()
                .map(|data| data / p_score_sum)
                .collect();
        }

        return (u_return, nu, refs, reads);
    }

    ///modifies the scores of u and nu with respect to max_score and min_score
    fn rescale_samscore(
        mut u: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        mut nu: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        max_score: f64,
        min_score: f64,
    ) -> (
        HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
    ) {
        let scaling_factor: f64;

        if min_score < 0.0 {
            scaling_factor = 100.0 / (max_score - min_score);
        } else {
            scaling_factor = 100.0 / max_score;
        }

        for k in u.clone().keys() {
            if min_score < 0.0 {
                u.get_mut(k).unwrap().1[0] = u.get(k).unwrap().1[0].clone() - min_score;
            }

            u.get_mut(k).unwrap().1[0] = f64::exp(u.get(k).unwrap().1[0] * scaling_factor);
            u.get_mut(k).unwrap().3 = u.get(k).unwrap().1[0];
        }

        for k in nu.clone().keys() {
            nu.get_mut(k).unwrap().3 = 0.0;

            for i in 0..nu.get(k).unwrap().1.len() {
                if min_score < 0.0 {
                    nu.get_mut(k).unwrap().1[i] = nu.get(k).unwrap().1[i] - min_score;
                }

                nu.get_mut(k).unwrap().1[i] = f64::exp(nu.get(k).unwrap().1[i] * scaling_factor);

                if nu.get(k).unwrap().1[i] > nu.get(k).unwrap().3 {
                    nu.get_mut(k).unwrap().3 = nu.get(k).unwrap().1[i];
                }
            }
        }
        return (u, nu);
    }
}

mod parse_sam {
    use std::{
        fmt::Debug,
        fs::File,
        io::{BufRead, BufReader},
    };

    ///Stores the desired fields of a .SAM record and the line itself as a String
    #[derive(Debug)]
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

    impl SamLine {
        ///Create a new Some(SamLine) object by consuming a String object
        ///
        ///Returns none if the provided String is to be ignored
        pub fn new(new_line: String) -> Option<SamLine> {
            if new_line.is_empty() || new_line.starts_with("#") || new_line.starts_with("@") {
                return None;
            }

            let fields = new_line.split("\t").collect::<Vec<&str>>();

            //extremely inefficient; should optimize later on
            let mut new_sam_line = SamLine {
                read_id: String::from(*(fields.get(0).expect("error parsing read_id"))),
                read_length: fields.get(9).expect("error parsing length field").len(),
                position: fields
                    .get(3)
                    .expect("error reading position field")
                    .parse::<u32>()
                    .expect("error parsing position as u32"),
                score: None,
                btws_flg: fields
                    .get(1)
                    .expect("error reading btws_flg field")
                    .parse::<u32>()
                    .expect("error parsing btws_flg as u32"),
                unmapped: ((fields.get(1).unwrap().parse::<u32>().unwrap()) & (4 as u32)
                    == (4 as u32)),
                ref_id: String::from(*(fields.get(2).expect("error parsing ref_id"))),
                sam_fields: fields
                    .into_iter()
                    .map(|data| String::from(data))
                    .collect::<Vec<String>>(),
                line: new_line,
            };

            new_sam_line.score = Some(find_sam_align_score(&mut new_sam_line));

            return Some(new_sam_line);
        }
    }

    fn find_sam_align_score(data: &SamLine) -> f64 {
        for field in data.sam_fields.clone() {
            if field.starts_with("AS:i:") {
                return (field[5..]
                    .parse::<i32>()
                    .expect("unable to parse field as i32 in find_sam_align_score")
                    as f64)
                    + (data.read_length as f64);
            }
        }

        panic!("unable to find sam alignment score!")
    }

    /// stores the result of parsing one line of a .SAM file\
    /// * Ok(T) => T is a SamLine object;\
    /// * Ignore and EOF are special flags;\
    /// * Err(String) => String indicates an error generating a SamLine object.
    pub enum ParseResult<T> {
        Ok(T),
        Ignore,
        EOF,
        Err(String),
    }

    pub fn parse_sam(
        sam_reader: &mut BufReader<File>,
        p_score_cutoff: Option<f64>,
    ) -> ParseResult<SamLine> {
        let p_score_cutoff: f64 = p_score_cutoff.unwrap_or(0.01);

        let mut buf: String = String::new();
        match sam_reader.read_line(&mut buf) {
            Ok(code) => {
                if code == 0 {
                    return ParseResult::EOF;
                } else {
                    match SamLine::new(buf) {
                        None => return ParseResult::Ignore,
                        Some(new_line) => {
                            if new_line
                                .score
                                .expect("error unwrapping newline.score in parseSAM")
                                > p_score_cutoff
                            {
                                return ParseResult::Ok(new_line);
                            } else {
                                return ParseResult::Ignore;
                            }
                        }
                    }
                }
            }
            Err(_) => {
                return ParseResult::Err(String::from(
                    "Error propagated in parseSAM from SamLine::new",
                ))
            }
        }
    }
}

pub fn compute_best_hit(
    u: &HashMap<i32, (i32, f64)>,
    nu: &HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
    refs: &Vec<String>,
    reads: &Vec<String>,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ref_count = refs.len();
    let mut best_hit_reads = vec![0.0; ref_count];
    let mut level_1_reads = vec![0.0; ref_count];
    let mut level_2_reads = vec![0.0; ref_count];

    for i in u.keys() {
        *(best_hit_reads
            .get_mut(u.get(i).unwrap().0 as usize)
            .unwrap()) += 1.0;
        *(level_1_reads.get_mut(u.get(i).unwrap().0 as usize).unwrap()) += 1.0;
    }

    for i in nu.keys() {
        let z = nu.get(i).unwrap();
        let ind = &z.0;
        let x_norm = &z.2;
        let best_ref = x_norm.iter().cloned().fold(-1. / 0. /* -inf */, f64::max);
        let mut num_best_ref = 0;

        for (j, _) in x_norm.iter().enumerate() {
            if *(x_norm.get(j).unwrap()) == best_ref {
                num_best_ref += 1;
            }
        }

        num_best_ref = match num_best_ref {
            0 => 1,
            _ => num_best_ref,
        };

        for (j, _) in x_norm.iter().enumerate() {
            if *(x_norm.get(j).unwrap()) == best_ref {
                *(best_hit_reads
                    .get_mut(*(ind.get(j).unwrap()) as usize)
                    .unwrap()) += 1.0 / num_best_ref as f64;

                if *(x_norm.get(j).unwrap()) >= 0.5 {
                    *(level_1_reads
                        .get_mut(*(ind.get(j).unwrap()) as usize)
                        .unwrap()) += 1.0;
                } else if *(x_norm.get(j).unwrap()) >= 0.01 {
                    *(level_2_reads
                        .get_mut(*(ind.get(j).unwrap()) as usize)
                        .unwrap()) += 1.0;
                }
            }
        }
    }

    let read_count = reads.len();

    let best_hit: Vec<f64> = best_hit_reads
        .iter()
        .map(|val| val.clone() / read_count as f64)
        .collect();
    let level1: Vec<f64> = level_1_reads
        .iter()
        .map(|val| val.clone() / read_count as f64)
        .collect();
    let level2: Vec<f64> = level_2_reads
        .iter()
        .map(|val| val.clone() / read_count as f64)
        .collect();

    return (best_hit_reads, best_hit, level1, level2);
}

pub fn em(
    u: &HashMap<i32, (i32, f64)>,
    mut nu: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
    genomes: &Vec<String>,
    max_iter: i32,
    epsilon: f64,
    pi_prior: f64,
    theta_prior: f64,
) -> (
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
) {
    let genome_count = genomes.len();
    let mut pi = vec![1.0 as f64 / genome_count as f64; genome_count];
    let mut init_pi = pi.clone();
    let mut theta = pi.clone();

    let mut pi_sum_0 = vec![0.0; genome_count];

    let u_weights: Vec<f64> = u.iter().map(|entry| (*(entry.1)).1).collect();
    let mut max_u_weights = 0.0;
    let mut u_total = 0.0;

    if !u_weights.is_empty() {
        max_u_weights = u_weights
            .iter()
            .cloned()
            .fold(-1. / 0. /* -inf */, f64::max);
        u_total = u_weights.iter().sum();
    }

    for i in u.keys() {
        pi_sum_0[u.get(i).unwrap().0 as usize] += u.get(i).unwrap().1;
    }

    let nu_weights: Vec<f64> = nu.iter().map(|entry| (*(entry.1)).3).collect();
    let mut max_nu_weights = 0.0;
    let mut nu_total = 0.0;

    if !nu_weights.is_empty() {
        max_nu_weights = nu_weights
            .iter()
            .cloned()
            .fold(-1. / 0. /* -inf */, f64::max);
        nu_total = nu_weights.iter().sum();
    }

    let prior_weight = f64::max(max_u_weights, max_nu_weights);
    let mut nu_length = nu.len();

    if nu_length == 0 {
        nu_length = 1;
    }

    //EM iterations
    for i in 0..max_iter {
        let pi_old = pi.clone();
        let mut theta_sum = vec![0.0; genome_count];

        //E step
        for j in nu.clone().keys() {
            let z = nu.get(j).unwrap().clone();

            //A set of any genome mapping with j
            let ind = &z.0;

            //Get relevant pis for the read
            let pi_temp: Vec<f64> = ind.iter().map(|val| pi[*val as usize].clone()).collect();

            //Get relevant thetas for the read
            let theta_temp: Vec<f64> = ind.iter().map(|val| theta[*val as usize].clone()).collect();

            //Calculate non-normalized xs
            let mut x_temp: Vec<f64> = Vec::new();

            for k in 0..ind.len() {
                x_temp.push(pi_temp[k] * theta_temp[k] * z.1[k]);
            }

            let x_sum: f64 = x_temp.iter().sum();

            //Avoid dividing by 0 at all times
            let x_norm: Vec<f64>;

            if x_sum == 0.0 {
                x_norm = vec![0.0; x_temp.len()];
            } else {
                x_norm = x_temp.iter().map(|val| val / x_sum).collect();
            }

            //Update x in nu
            nu.get_mut(j).unwrap().2 = x_norm.clone();

            for (k, _) in ind.iter().enumerate() {
                theta_sum[ind[k] as usize] += x_norm[k] * nu.get(j).unwrap().3;
            }
        }

        //M step
        let pi_sum: Vec<f64> = theta_sum
            .iter()
            .enumerate()
            .map(|(idx, _)| theta_sum[idx] + pi_sum_0[idx])
            .collect();
        let pip = pi_prior * prior_weight;

        //update pi
        pi = pi_sum
            .iter()
            .map(|val| ((*val as f64) + pip) / (u_total + nu_total + (pip * pi_sum.len() as f64)))
            .collect();

        if i == 0 {
            init_pi = pi.clone();
        }

        let theta_p = theta_prior * prior_weight;

        let mut nu_total_div = nu_total;

        if nu_total_div == 0 as f64 {
            nu_total_div = 1 as f64;
        }

        theta = theta_sum
            .iter()
            .map(|val| (*val + theta_p) / (nu_total_div + (theta_p * theta_sum.len() as f64)))
            .collect();

        let mut cutoff = 0.0;

        for (k, _) in pi.iter().enumerate() {
            cutoff += (pi_old[k] - pi[k]).abs();
        }

        if cutoff <= epsilon || nu_length == 1 {
            break;
        }
    }

    return (init_pi, pi, theta, nu);
}

mod rewrite_align {
    use crate::parse_sam::*;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;
    use std::io::LineWriter;
    use std::io::Write;

    pub fn rewrite_align(
        u: &HashMap<i32, (i32, f64)>,
        nu: &HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        sam_path: &str,
        p_score_cutoff: &f64,
        path: &str,
    ) {
        let mut read_id_dict: HashMap<String, i32> = HashMap::new();
        let mut ref_id_dict: HashMap<String, i32> = HashMap::new();

        let mut genomes: Vec<String> = Vec::new();
        let mut read: Vec<String> = Vec::new();

        let mut ref_count = 0;
        let mut read_count = 0;

        let old_file = File::open(sam_path).expect("Invalid file");
        let mut sam_reader = BufReader::new(old_file);
        let new_file = File::create(path).expect("unable to create file");
        let mut sam_writer = LineWriter::new(new_file);

        //for line in parseSam
        loop {
            let sam_line = parse_sam(&mut sam_reader, Some(*p_score_cutoff));

            let sam_line = match sam_line {
                ParseResult::Ok(line) => line,
                ParseResult::Ignore => continue,
                ParseResult::EOF => break,
                ParseResult::Err(_) => {
                    panic!("unable to read old_file in rewrite_align::rewrite_align")
                }
            };

            let mut ref_index = ref_id_dict.get(&sam_line.ref_id).unwrap_or(&-1).clone();

            if ref_index == -1 {
                ref_index = ref_count.clone();
                ref_id_dict.insert(sam_line.ref_id.clone(), ref_index);
                genomes.push(sam_line.ref_id);
                ref_count += 1;
            }

            let mut read_index = read_id_dict.get(&sam_line.read_id).unwrap_or(&-1).clone();

            if read_index == -1 {
                // hold on to this new read
                // first, wrap previous read profile and see if any previous read has a
                // same profile with that!

                read_index = read_count.clone();
                read_id_dict.insert(sam_line.read_id.clone(), read_index);
                read.push(sam_line.read_id);
                read_count += 1;

                if u.contains_key(&read_index) {
                    sam_writer
                        .write(sam_line.line.as_bytes())
                        .expect("unable to write to new_file in rewrite_align::rewrite_align");
                    continue;
                }
            }

            if nu.contains_key(&read_index) {
                if find_updated_score(&nu, read_index, ref_index) < *p_score_cutoff {
                    continue;
                }
                sam_writer
                    .write(sam_line.line.as_bytes())
                    .expect("unable to write to new_file in rewrite_align::rewrite_align");
            }
        }
    }

    fn find_updated_score(
        nu: &HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        read_index: i32,
        ref_index: i32,
    ) -> f64 {
        let v1 = match nu.get(&read_index) {
            Some(val) => val,
            None => return 0.0,
        };

        let mut idx: usize = 0;

        for (i, el) in v1.0.iter().enumerate() {
            if *el == ref_index {
                idx = i;
                break;
            }
        }

        return nu.get(&read_index).unwrap().2.get(idx).unwrap().clone();
    }
}

///tests and whatnot
#[cfg(test)]
mod tests {

    #![allow(unused)]

    use crate::build_matrix::*;
    use crate::rewrite_align::*;
    use crate::*;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    use std::io::Read;

    extern crate yaml_rust;
    use yaml_rust::{YamlEmitter, YamlLoader};

    #[test]
    fn test_rewrite_align() {
        let (u, nu, refs, reads) = build_matrix("TestFiles/test_al.sam", None);
        let (init_pi, pi, theta, nu) = em(&u, nu, &refs, 5, 1e-7, 0.0, 0.0);
        rewrite_align(
            &u,
            &nu,
            "TestFiles/test_al.sam",
            &0.01,
            "TestFiles/rewrite.sam",
        );

        let mut new_file =
            BufReader::new(File::open("TestFiles/rewrite.sam").expect("Invalid file"));
        let mut test_file = BufReader::new(
            File::open("tests/test_pathoscope/test_rewrite_align.txt").expect("Invalid file"),
        );

        let mut new_line = String::new();
        let mut test_line = String::new();

        //compare output
        loop {
            new_line.clear();
            test_line.clear();

            match new_file.read_line(&mut new_line) {
                Ok(val) => {
                    if val == 0 {
                        return;
                    }
                }

                Err(_) => panic!("tests::test_rewrite_align: `error reading rewrite.sam`"),
            }

            match test_file.read_line(&mut test_line) {
                Ok(val) => {
                    if val == 0 {
                        return;
                    }
                }

                Err(_) => {
                    panic!("tests::test_rewrite_align `error reading test_rewrite_align.txt`")
                }
            }

            let new_line_fields: Vec<&str> = new_line.split('\t').collect();
            let test_line_fields: Vec<&str> = test_line.split('\t').collect();

            //compare fields of desired output to actual output
            for i in 0..new_line_fields.len() {
                assert!(
                    new_line_fields[i] == test_line_fields[i],
                    "tests::test_rewrite_align: `output of rewrite_align does not match expected`"
                );
            }
        }
    }

    #[test]
    fn test_em() {
        let (u, nu, refs, reads) = build_matrix("TestFiles/test_al.sam", None);
        let (init_pi, pi, theta, nu) = em(&u, nu, &refs, 5, 1e-6, 0.0, 0.0);

        let mut test_file = File::open("tests/test_pathoscope/test_em_5_1e_06_0_0_.yml")
            .expect("tests::test_build_matrix: `unable to open test file`");
        let mut test_string = String::new();
        test_string.clear();
        test_file
            .read_to_string(&mut test_string)
            .expect("tests::test_build_matrix: unable to read test file");
        let test_matrix = YamlLoader::load_from_str(&test_string)
            .expect("tests::test_build_matrix: unable to parse test file as .yml")[0]
            .clone();

        //compare nu
        for (key, value) in nu {
            for i in 0..value.0.len() {
                assert!(value.0[i] as i64 == test_matrix[3][key as usize][0][i].as_i64().unwrap());
            }
        }

        let mut init_pi_count = 0;

        //compare init_pi
        for entry in &init_pi {
            //vector is not sorted; check every index and break if found
            for i in 0..init_pi.len() {
                if ((*entry) - test_matrix[0][i].as_f64().unwrap()) <= 0.000000000000001 {
                    init_pi_count += 1;
                    break;
                } else {
                    continue;
                }
            }
        }
        if init_pi_count != test_matrix[0].as_vec().unwrap().len() {
            panic!();
        }

        let mut pi_count = 0;

        //compare pi
        for entry in &pi {
            //vector is not sorted; check every index and break if found
            for i in 0..pi.len() {
                if ((*entry) - test_matrix[1][i].as_f64().unwrap()) <= 0.000000000000001 {
                    pi_count += 1;
                    break;
                } else {
                    continue;
                }
            }
        }
        if pi_count != test_matrix[1].as_vec().unwrap().len() {
            panic!();
        }

        let mut theta_count = 0;

        //compare pi
        for entry in &theta {
            //vector is not sorted; check every index and break if found
            for i in 0..theta.len() {
                if ((*entry) - test_matrix[2][i].as_f64().unwrap()) <= 0.000000000000001 {
                    theta_count += 1;
                    break;
                } else {
                    continue;
                }
            }
        }
        if theta_count != test_matrix[2].as_vec().unwrap().len() {
            panic!();
        }
    }

    #[test]
    fn test_best_hit() {
        let (u, nu, refs, reads) = build_matrix("TestFiles/test_al.sam", None);
        let (best_hit_reads, best_hit, level1, level2) = compute_best_hit(&u, &nu, &refs, &reads);

        let mut test_file = File::open("tests/test_pathoscope/test_compute_best_hit.yml")
            .expect("tests::test_build_matrix: `unable to open test file`");
        let mut test_string = String::new();
        test_string.clear();
        test_file
            .read_to_string(&mut test_string)
            .expect("tests::test_build_matrix: unable to read test file");
        let test_matrix = YamlLoader::load_from_str(&test_string)
            .expect("tests::test_build_matrix: unable to parse test file as .yml")[0]
            .clone();

        //compare best_hit_reads
        for i in 0..best_hit_reads.len() {
            assert!(best_hit_reads[i] == test_matrix[0][i].as_f64().unwrap())
        }

        //compare best_hit
        for i in 0..best_hit.len() {
            assert!(best_hit[i] == test_matrix[1][i].as_f64().unwrap());
        }

        //compare level1
        for i in 0..level1.len() {
            assert!(level1[i] == test_matrix[2][i].as_f64().unwrap());
        }

        //compare level2
        for i in 0..level2.len() {
            assert!(level2[i] == test_matrix[3][i].as_f64().unwrap());
        }
    }

    #[test]
    fn test_build_matrix() {
        let (u, nu, refs, reads) = build_matrix("TestFiles/test_al.sam", None);

        let mut test_file = File::open("tests/test_pathoscope/test_build_matrix.yml")
            .expect("tests::test_build_matrix: `unable to open test file`");
        let mut test_string = String::new();
        test_string.clear();
        test_file
            .read_to_string(&mut test_string)
            .expect("tests::test_build_matrix: unable to read test file");
        let test_matrix = YamlLoader::load_from_str(&test_string)
            .expect("tests::test_build_matrix: unable to parse test file as .yml")[0]
            .clone();

        //compare u
        for (key, value) in u {
            assert!(value.0 as i64 == test_matrix[0][key as usize][0].as_i64().unwrap());
            assert!(value.1 == test_matrix[0][key as usize][1].as_f64().unwrap());
        }

        //compare nu
        for (key, value) in nu {
            for i in 0..value.0.len() {
                assert!(value.0[i] as i64 == test_matrix[1][key as usize][0][i].as_i64().unwrap());
            }

            for i in 0..value.1.len() {
                assert!(value.1[i] == test_matrix[1][key as usize][1][i].as_f64().unwrap());
            }

            for i in 0..value.2.len() {
                assert!(value.2[i] == test_matrix[1][key as usize][2][i].as_f64().unwrap());
            }

            assert!(value.3 == test_matrix[1][key as usize][3].as_f64().unwrap());
        }

        let mut ref_count = 0;

        //compare refs
        for entry in &refs {
            //vector is not sorted; check every index and break if found
            for i in 0..refs.len() {
                if (*entry).eq(test_matrix[2][i].as_str().unwrap()) {
                    ref_count += 1;
                    break;
                } else {
                    continue;
                }
            }
        }
        if ref_count != test_matrix[2].as_vec().unwrap().len() {
            panic!();
        }

        let mut read_count = 0;

        //compare reads
        for entry in &reads {
            //vector is not sorted; check every index and break if found
            for i in 0..reads.len() {
                if (*entry).eq(test_matrix[3][i].as_str().unwrap()) {
                    read_count += 1;
                    break;
                } else {
                    continue;
                }
            }
        }
        if read_count != test_matrix[3].as_vec().unwrap().len() {
            panic!();
        }
    }
}
