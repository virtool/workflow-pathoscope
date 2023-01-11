#![allow(nonstandard_style)]
#![allow(unused)]


///Wrapper for the buildMatrix function and adjacent code
mod buildMatrix
{
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;
    use super::parseSam::*;

    ///produces the SAM matrix given a path to a .SAM file
    pub fn buildMatrix(SAMPath: &str, pScoreCutoff: Option<f64>) -> (
                                                                    HashMap<i32, (i32, f64)>, 
                                                                    HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
                                                                    Vec<String>,
                                                                    Vec<String>
                                                                )
    {
        let mut u: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)> = HashMap::new();
        let mut nu: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)> = HashMap::new();

        let mut hReadId: HashMap<String, i32> = HashMap::new();
        let mut hRefId: HashMap<String, i32> = HashMap::new();

        let mut refs: Vec<String> = Vec::new();
        let mut reads: Vec<String> = Vec::new();

        let mut refCount: i32 = 0;
        let mut readCount: i32 = 0;

        let mut maxScore: f64 = 0.0;
        let mut minScore: f64 = 0.0;


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
                        panic!();
                    }
                    parseResult::Ignore => continue
            }

            minScore = newLine.score.unwrap().min(minScore); // might need to completely redo this..
            maxScore = newLine.score.unwrap().max(maxScore);

            let mut refIndex = (hRefId.get(&newLine.ref_id).unwrap_or(&-1)).clone();

            if refIndex == -1
            {
                refIndex = refCount;
                hRefId.insert(newLine.ref_id.clone(), refIndex);
                refs.push(newLine.ref_id.clone());
                refCount += 1;
            }

            let mut readIndex = (hReadId.get(&newLine.read_id).unwrap_or(&-1)).clone();

            if readIndex == -1
            {
                readIndex = readCount;
                hReadId.insert(newLine.read_id.clone(), readIndex);
                reads.push(newLine.read_id.clone());
                readCount += 1;

                u.insert(readIndex, (
                    vec![refIndex],
                    vec![newLine.score.clone().unwrap()],
                    vec![newLine.score.clone().unwrap() as f64],
                    newLine.score.clone().unwrap()
                ));
            }
            else
            {
                if u.contains_key(&readIndex)
                {
                    if u.get(&readIndex).unwrap().0.contains(&refIndex)
                    {
                        continue;
                    }
                    nu.insert(readIndex, u.get(&readIndex).unwrap().clone());
                    u.remove(&readIndex);
                }

                if nu.get(&readIndex).unwrap().0.contains(&refIndex)
                {
                    continue;
                }

                nu.get_mut(&readIndex).unwrap().0.push(refIndex);
                nu.get_mut(&readIndex).unwrap().1.push(newLine.score.unwrap());

                if newLine.score.unwrap() > nu.get(&readIndex).unwrap().3
                {
                    nu.get_mut(&readIndex).unwrap().3 = newLine.score.unwrap();
                }
            }
        
        }

        let (mut u, mut nu) = rescaleSamscore(u, nu, maxScore, minScore);
        
        let mut uPtr = &mut u as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;
        let mut nuPtr = &mut nu as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;

        let mut uReturn: HashMap<i32, (i32, f64)> = HashMap::new();
        let mut uReturnPtr = &mut uReturn as *mut HashMap<i32, (i32, f64)>;

        unsafe
        {
        for k in (*uPtr).keys()
        {
            (*uReturnPtr).insert(k.clone(), (
                (*uPtr).get(k).unwrap().0.get(0).unwrap().clone(),
                (*uPtr).get(k).unwrap().1.get(0).unwrap().clone()
            ));
        }

        for k in (*nuPtr).keys()
        {
            let pScoreSum = (*nuPtr).get(k).unwrap().1.iter().sum::<f64>();

                (*nuPtr).get_mut(k).unwrap().2 = (*nuPtr).get(k).unwrap().1.iter().map(|data| {data / pScoreSum}).collect();
        }
        }

        return (uReturn, nu, refs, reads);
        
    }

    ///modifies the scores of u and nu with respect to maxScore and minScore
    fn rescaleSamscore  ( 
                            mut u: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
                            mut nu: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
                            maxScore: f64,
                            minScore: f64
                        ) ->
                        (
                            HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
                            HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>
                        )
    {
        let mut scalingFactor: f64;

        //Need to access the data 3 places at once so this is a quick solution;
        //will possibly replace with smart pointers later on
        let mut uPtr = &mut u as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;
        let mut nuPtr = &mut nu as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;

        if minScore < 0.0
        {
            scalingFactor = 100.0 / (maxScore - minScore);
        }
        else
        {
            scalingFactor = 100.0 / maxScore;
        }

        //Derefing raw pointers lmao
        unsafe
        {
        for k in (*uPtr).keys()
        {
            if minScore < 0.0
            {
                (*uPtr).get_mut(k).unwrap().1[0] = (*uPtr).get(k).unwrap().1[0].clone() - minScore;
            }

            (*uPtr).get_mut(k).unwrap().1[0] = f64::exp((*uPtr).get(k).unwrap().1[0] * scalingFactor);
            (*uPtr).get_mut(k).unwrap().3 = u.get(k).unwrap().1[0];
        }

        for k in (*nuPtr).keys()
        {
            (*nuPtr).get_mut(k).unwrap().3 = 0.0; //idk what this does; im just copy pasting it from python tbh

            for i in (0..(*nuPtr).get(k).unwrap().1.len())
            {
                if minScore < 0.0
                {
                    (*nuPtr).get_mut(k).unwrap().1[i] = ((*nuPtr).get(k).unwrap().1[i] - minScore);
                }

                (*nuPtr).get_mut(k).unwrap().1[i] = f64::exp((*nuPtr).get(k).unwrap().1[i] * scalingFactor);

                if (*nuPtr).get(k).unwrap().1[i] > (*nuPtr).get(k).unwrap().3
                {
                    (*nuPtr).get_mut(k).unwrap().3 = (*nuPtr).get(k).unwrap().1[i];
                }
            }
        }
        }
        return (u, nu);
    }
}


///Wrapper for the parseSAM function and adjacent code
mod parseSam
{
    use std::{io::{BufReader, BufRead}, fs::File, fmt::Debug, rc::Rc};

    #[derive(Debug)]
    pub struct SamLine
    {
        pub read_id: String,
        pub read_length: usize,
        pub position: u32,
        pub score: Option<f64>,
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


    fn find_sam_align_score(data: &SamLine) -> f64
    {
        let readLength: f64 = data.read_length as f64;

        for field in data.SAMfields.clone()
        {
            if field.starts_with("AS:i:")
            {
                return (field[5..].parse::<i32>().expect("unable to parse field as i32 in find_sam_align_score") as f64) + (data.read_length as f64);
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
    pub fn parseSAM(SAMReader: &mut BufReader<File>, pScoreCutoff: Option<f64>) -> parseResult<SamLine>
    {
        let pScoreCutoff: f64 = pScoreCutoff.unwrap_or(0.01);

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
                    match SamLine::new(buf)
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


///tests and whatnot
#[cfg(test)]
mod tests
{
    use std::io::BufReader;
    use std::fs::File;

    use crate::parseSam::*;
    use crate::buildMatrix::*;

    ///tests the buildMatrix function
    #[test]
    fn testBuildMatrix()
    {
        let (u, nu, refs, reads) = buildMatrix("TestFiles/test_al.sam", None);
        
        //testing u length and values
        assert!(u.len() == 19174);
       
        assert!(u.get(&0).unwrap().0 == 0);
        assert!(u.get(&0).unwrap().1 == 1.3892652283160566e+43);


        //testing nu length and values
        assert!(nu.len() == 1102);

        assert!(nu.get(&10).unwrap().0[0] == 1);
        assert!(nu.get(&10).unwrap().0[1] == 3);
        assert!(nu.get(&10).unwrap().0.len() == 2);

        assert!(nu.get(&10).unwrap().1[0] == 9.539599502409837e+36);
        assert!(nu.get(&10).unwrap().1[1] == 9.394830991271991e+34);
        assert!(nu.get(&10).unwrap().1.len() == 2);

        assert!(nu.get(&10).unwrap().2[0] == 0.9902477974114013);
        assert!(nu.get(&10).unwrap().2[1] == 0.009752202588598546);
        assert!(nu.get(&10).unwrap().2.len() == 2);
        
        assert!(nu.get(&10).unwrap().3 == 9.539599502409837e+36);
        

        //testing refs length and values
        assert!(refs.len() == 40);

        assert!(refs.get(0).unwrap().eq("NC_016509"));
        assert!(refs.get(1).unwrap().eq("NC_003615"));
        assert!(refs.get(2).unwrap().eq("NC_007448"));


        //testing reads length and values
        assert!(reads.len() == 20276);

        assert!(reads.get(0).unwrap().eq("HWI-ST1410:82:C2VAGACXX:7:1101:20066:1892"));
        assert!(reads.get(1).unwrap().eq("HWI-ST1410:82:C2VAGACXX:7:1101:11037:2144"));
        assert!(reads.get(2).unwrap().eq("HWI-ST1410:82:C2VAGACXX:7:1101:14679:2757"));
        
    }


    ///tests the parseSAM function
    #[test]
    fn testParseSAM()
    {
        println!("Testing parseSAM");

        let SAMFile = File::open("TestFiles/test_al.sam").expect("Invalid file");
        let mut SAMReader = BufReader::new(SAMFile);
        loop
        {
            let newLine = parseSAM(&mut SAMReader, None);
            match newLine
            {
                parseResult::Ok(data) => continue,
                parseResult::EOF => return,
                parseResult::Err(msg) => return,
                parseResult::Ignore => continue
            }
        }
    }
}