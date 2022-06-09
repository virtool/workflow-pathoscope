use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::Write;
use std::io::{self, BufRead};
use std::path::Path;

fn main() {
    let args: Vec<String> = env::args().collect();

    let isolate_sam_path = &args[1];
    let subtraction_sam_path = &args[2];
    let output_sam_path = &args[3];

    let subtraction_scores = parse_subtraction_sam(subtraction_sam_path);

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

// Check if the passed read_id should be eliminated if its isolate score is
// higher than the subtraction score.
fn check_should_eliminate(
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
fn find_sam_align_score(fields: &Vec<&str>) -> f32 {
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

fn parse_subtraction_sam(path: &String) -> HashMap<String, f32> {
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

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
