#![allow(nonstandard_style)]

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
    ///Returns none if the provided String is empty
    pub fn new(newLine: String) -> Option<SamLine>
    {
        if newLine.is_empty()
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
        //parse line into fields
        newSamLine.fields = Some(newSamLine.line.split("\t").map(|s| {String::from(s)}).collect());


/*To do here..:
    update the constructor of the samline class so its not dumb like it is rn
    set the rest of the class attributes
*/







        return Some(newSamLine);
    }    

}





///parses one line of a .SAM file, taking reference to a BufReader<File>
/// 
///returns Some(SamLine) if data was read or None if none was read
fn parseSam(SAMReader: &mut BufReader<File>) -> Option<SamLine>
{
    let mut buf: String = String::new();
    match SAMReader.read_line(&mut buf)
    {
        Ok(_) => return SamLine::new(String::from(buf.trim())),
        Err(_) => return None
    }   
}


#[cfg(test)]
mod tests {
    use super::*;

    ///tests the parseSam function
    #[test]
    fn testParseSam() {
        let SAMFile = File::open("TestFiles/to_isolates.sam").expect("Invalid file");
        let mut SAMReader = BufReader::new(SAMFile);
        loop
        {
            let newLine = super::parseSam(&mut SAMReader);
            match newLine
            {
                Some(data) => println!("{:?}", data.fields),
                None => return
            }
        }
    }
}