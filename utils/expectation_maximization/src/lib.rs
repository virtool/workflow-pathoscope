#![allow(nonstandard_style)]
#![allow(unused)]

use std::{io::{BufReader, BufRead}, fs::File};

struct SamLine
{
    line: String,

    fields: Option<Vec<String>>,

    read_id: Option<String>,
    read_length: Option<u32>,
    position: Option<u32>,
    score: Option<f32>,
    btwsFlg: Option<u32>,
    unmapped: Option<bool>,
    ref_id: Option<String>,
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

        let mut newSamLine = SamLine
        { 
            line: newLine,

            fields: None,

            read_id: None, 
            read_length: None, 
            position: None, 
            score: None, 
            btwsFlg: None, 
            unmapped: None, 
            ref_id: None
        };

        newSamLine.fields = Some(newSamLine.line.split("\t").map(|s| {String::from(s)}).collect());

        newSamLine.read_id = Some(String::from(newSamLine.fields.as_ref().unwrap()[0].as_str()));
        newSamLine.read_length = Some(newSamLine.fields.as_ref().unwrap()[9].len() as u32);
        newSamLine.position = Some(newSamLine.fields.as_ref().unwrap()[3].parse::<u32>().unwrap());
        newSamLine.btwsFlg = Some(newSamLine.fields.as_ref().unwrap()[1].parse::<u32>().unwrap());
        newSamLine.unmapped = Some(newSamLine.btwsFlg.as_ref().unwrap() & 4 == 4);
        newSamLine.ref_id = Some(String::from(newSamLine.fields.as_ref().unwrap()[2].as_str()));

        newSamLine.score = None; //implement later!


        return Some(newSamLine);
    }    
}


fn find_sam_align_score(fields: Option<Vec<String>>)
{
    unimplemented!()
}


enum parseResult<T>
{
    Ok(T),
    Ignore,
    EOF,
    Err(String)
}

///parses one line of a .SAM file, taking reference to a BufReader<File>
/// 
///returns Some(SamLine) if data was read or None if none was read
fn parseSam(SAMReader: &mut BufReader<File>) -> parseResult<SamLine>
{
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
                    Some(newLine) => return parseResult::Ok(newLine)
                }
            }
            
        }
        Err(_) => return parseResult::Err(String::from("Error propagated in parseSam from SamLine::new"))
    }   
}


#[cfg(test)]
mod tests {
    use super::*;

    ///tests the parseSam function
    #[test]
    fn testParseSam() {
        let SAMFile = File::open("TestFiles/test_al.sam").expect("Invalid file");
        let mut SAMReader = BufReader::new(SAMFile);
        loop
        {
            let newLine = super::parseSam(&mut SAMReader);
            match newLine
            {
                parseResult::Ok(data) => println!("{:?}", data.fields),
                parseResult::EOF => return,
                parseResult::Err(msg) => return,
                parseResult::Ignore => continue
            }
        }
    }
}