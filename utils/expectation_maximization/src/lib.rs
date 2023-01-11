#![allow(nonstandard_style)]
#![allow(unused)]

use std::{path::Path, collections::{hash_map, HashMap}, fs::File, io::BufReader, cmp::min};
use crate::SAMParser::*;

use SAMParser::parseSAM;


fn buildMatrix(SAMPath: &str, pScoreCutoff: Option<f32>)
{
    let mut u: HashMap<i32, (Vec<i32>, Vec<f32>, Vec<f32>, f32)> = HashMap::new();
    let mut nu: HashMap<i32, ((i32, i32), (f32, f32), f32, f32)> = HashMap::new();

    let mut hReadId: HashMap<String, i32> = HashMap::new();
    let mut hRefId: HashMap<String, i32> = HashMap::new();

    let mut refs: Vec<String> = Vec::new();
    let mut reads: Vec<String> = Vec::new();

    let mut refCount: i32 = 0;
    let mut readCount: i32 = 0;

    let mut maxScore: f32 = 0.0;
    let mut minScore: f32 = 0.0;


    let SAMFile = File::open("TestFiles/test_al.sam").expect("Invalid file");
    let mut SAMReader = BufReader::new(SAMFile);
    
    loop
    {
        let newLine: SamLine;
        
        match parseSAM(&mut SAMReader, pScoreCutoff)
        {
                parseResult::Ok(data) => 
                {
                    if data.score < pScoreCutoff
                    {
                        continue;
                    }
                    else
                    {
                        newLine = data;
                    }
                }
                parseResult::EOF => break,
                parseResult::Err(msg) => 
                {
                    println!("{}",msg);
                    return;
                }
                parseResult::Ignore => continue
        }

        minScore = newLine.score.unwrap().min(minScore);
        maxScore = newLine.score.unwrap().max(maxScore);

        let mut refIndex = (hRefId.get(&newLine.ref_id).unwrap_or(&-1)).clone();

        if refIndex == -1
        {
            refIndex = refCount;
            hRefId.insert(newLine.ref_id.clone(), refIndex);
            refs.push(newLine.ref_id.clone());
            refCount += 1;
        }
    }

    println!("boop!");
}



///Wrapper for the parseSAM function, which serves as the entry point for the module.
mod SAMParser
{
    use std::{io::{BufReader, BufRead}, fs::File, fmt::Debug};

    #[derive(Debug)]
    pub struct SamLine
    {
        pub read_id: String,
        pub read_length: usize,
        pub position: u32,
        pub score: Option<f32>,
        pub btwsFlg: u32,
        pub unmapped: bool,
        pub ref_id:String,
        pub SAMfields: Vec<String>
    }


    impl SamLine
    {
        ///Create a new SamLine object by consuming a String object
        ///
        ///Returns none if the provided String is to be ignored
        pub fn new(newLine: String) -> Option<SamLine>
        {
            if newLine.is_empty() || newLine.starts_with("#") || newLine.starts_with("@")
            {
                return None;
            }

            let fields = newLine.split("\t").collect::<Vec<&str>>();

            let mut newSamLine = SamLine
            {
                read_id: String::from(*(fields.get(0).expect("error parsing read_id"))), 
                read_length: fields.get(9).expect("error parsing length field").len(),
                position: fields.get(3).expect("error reading position field").parse::<u32>().expect("error parsing position as u32"),
                score: None, 
                btwsFlg: fields.get(1).expect("error reading btwsFlg field").parse::<u32>().expect("error parsing btwsFlg as u32"), 
                unmapped: ((fields.get(1).unwrap().parse::<u32>().unwrap()) & (4 as u32) == (4 as u32)),
                ref_id: String::from(*(fields.get(2).expect("error parsing ref_id"))),
                SAMfields: fields.into_iter().map(|data| {String::from(data)}).collect::<Vec<String>>()
            };

            newSamLine.score = Some(find_sam_align_score(&mut newSamLine));

            return Some(newSamLine);
        }    
    }


    fn find_sam_align_score(data: &SamLine) -> f32
    {
        let readLength: f32 = data.read_length as f32;

        for field in data.SAMfields.clone()
        {
            if field.starts_with("AS:i:")
            {
                return (field[5..].parse::<i32>().expect("unable to parse field as i32 in find_sam_align_score") as f32) + (data.read_length as f32);
            }
        }

        panic!("unable to find sam alignment score!")


    }


    pub enum parseResult<T>
    {
        Ok(T),
        Ignore,
        EOF,
        Err(String)
    }

    ///parses one line of a .SAM file, taking reference to a BufReader<File>
    /// 
    ///returns Some(SamLine) if data was read or None if none was read
    pub fn parseSAM(SAMReader: &mut BufReader<File>, pScoreCutoff: Option<f32>) -> parseResult<SamLine>
    {
        let pScoreCutoff: f32 = pScoreCutoff.unwrap_or(0.01);

        let mut buf: String = String::new();
        match SAMReader.read_line(&mut buf)
        {
            Ok(code) => 
            {
                if code == 0
                {
                    return parseResult::EOF
                }
                else
                {
                    match SamLine::new(String::from(buf.trim()))
                    {
                        None => return parseResult::Ignore,
                        Some(newLine) =>
                        {
                            if newLine.score.expect("error unwrapping newline.score in parseSAM") > pScoreCutoff
                            {
                                return parseResult::Ok(newLine); 
                            }
                            else
                            {
                                return parseResult::Ignore;
                            }
                            
                        }
                    }
                }
                
            }
            Err(_) => return parseResult::Err(String::from("Error propagated in parseSAM from SamLine::new"))
        }   
    }
}


#[cfg(test)]
mod tests {
    use crate::SAMParser::*;
    use std::io::BufReader;
    use std::fs::File;
    use super::*;


    #[test]
    fn testBuildMatrix()
    {
        buildMatrix("TestFiles/test_al.sam", None)
    }



    ///tests the parseSAM function
    #[test]
    fn testParseSAM()
    {
        let SAMFile = File::open("TestFiles/test_al.sam").expect("Invalid file");
        let mut SAMReader = BufReader::new(SAMFile);
        loop
        {
            let newLine = parseSAM(&mut SAMReader, None);
            match newLine
            {
                parseResult::Ok(data) => println!("{:?}", data),
                parseResult::EOF => return,
                parseResult::Err(msg) => return,
                parseResult::Ignore => continue
            }
        }
    }
}