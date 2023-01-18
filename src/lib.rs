#![allow(nonstandard_style)]

///Wrapper for the buildMatrix function and adjacent code
mod BuildMatrix
{
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;
    use super::ParseSam::*;

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


        let SAMFile = File::open(SAMPath).expect("Invalid file");
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
        
        let uPtr = &mut u as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;
        let nuPtr = &mut nu as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;

        let mut uReturn: HashMap<i32, (i32, f64)> = HashMap::new();
        let uReturnPtr = &mut uReturn as *mut HashMap<i32, (i32, f64)>;

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
        let scalingFactor: f64;

        //Need to access the data 3 places at once so this is a quick solution;
        //will possibly replace with smart pointers later on
        let uPtr = &mut u as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;
        let nuPtr = &mut nu as *mut HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>;

        if minScore < 0.0
        {
            scalingFactor = 100.0 / (maxScore - minScore);
        }
        else
        {
            scalingFactor = 100.0 / maxScore;
        }

        //Derefing raw pointers
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
            (*nuPtr).get_mut(k).unwrap().3 = 0.0;

            for i in 0..(*nuPtr).get(k).unwrap().1.len()
            {
                if minScore < 0.0
                {
                    (*nuPtr).get_mut(k).unwrap().1[i] = (*nuPtr).get(k).unwrap().1[i] - minScore;
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
mod ParseSam
{
    use std::{io::{BufReader, BufRead}, fs::File, fmt::Debug};

    ///Stores the desired fields of a .SAM record and the line itself as a String
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
        pub SAMfields: Vec<String>,
        pub line: String
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

            //extremely inefficient; should optimize later on
            let mut newSamLine = SamLine
            {
                read_id: String::from(*(fields.get(0).expect("error parsing read_id"))), 
                read_length: fields.get(9).expect("error parsing length field").len(),
                position: fields.get(3).expect("error reading position field").parse::<u32>().expect("error parsing position as u32"),
                score: None, 
                btwsFlg: fields.get(1).expect("error reading btwsFlg field").parse::<u32>().expect("error parsing btwsFlg as u32"), 
                unmapped: ((fields.get(1).unwrap().parse::<u32>().unwrap()) & (4 as u32) == (4 as u32)),
                ref_id: String::from(*(fields.get(2).expect("error parsing ref_id"))),
                SAMfields: fields.into_iter().map(|data| {String::from(data)}).collect::<Vec<String>>(),
                line: newLine
            };

            newSamLine.score = Some(find_sam_align_score(&mut newSamLine));

            return Some(newSamLine);
        }    
    }


    fn find_sam_align_score(data: &SamLine) -> f64
    {

        for field in data.SAMfields.clone()
        {
            if field.starts_with("AS:i:")
            {
                return (field[5..].parse::<i32>().expect("unable to parse field as i32 in find_sam_align_score") as f64) + (data.read_length as f64);
            }
        }

        panic!("unable to find sam alignment score!")


    }

    /// stores the result of parsing one line of a .SAM file\
    /// * Ok(T) => T is a SamLine object;\
    /// * Ignore and EOF are special flags;\
    /// * Err(String) => String indicates an error generating a SamLine object.
    pub enum parseResult<T>
    {
        Ok(T),
        Ignore,
        EOF,
        Err(String)
    }


    ///parses one line of a .SAM file, taking reference to a BufReader<File>\
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

///wrapper for the computeBestHit function
mod BestHit
{
    use std::collections::HashMap;

    pub fn computeBestHit
    (  
        u: &HashMap<i32, (i32, f64)>,
        nu: &HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        refs: &Vec<String>,
        reads: &Vec<String>
    )
    ->
    (
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>
    )
    {
        let refCount = refs.len();
        let mut bestHitReads = vec![0.0; refCount];
        let mut level1Reads = vec![0.0; refCount];
        let mut level2Reads = vec![0.0; refCount];

        for i in u.keys()
        {
            *(bestHitReads.get_mut(u.get(i).unwrap().0 as usize).unwrap()) += 1.0;
            *(level1Reads.get_mut(u.get(i).unwrap().0 as usize).unwrap()) += 1.0;
        }

        for i in nu.keys()
        {
            let z = nu.get(i).unwrap();
            let ind = &z.0;
            let xNorm = &z.2;
            let bestRef = xNorm.iter().cloned().fold(-1./0. /* -inf */, f64::max);
            let mut numBestRef = 0;

            for (j, _) in xNorm.iter().enumerate()
            {
                if *(xNorm.get(j).unwrap()) == bestRef
                {
                    numBestRef += 1;
                }
            }

            numBestRef = match numBestRef
            {
                0 => 1,
                _ => numBestRef
            };

            
            for (j, _) in xNorm.iter().enumerate()
            {
                if *(xNorm.get(j).unwrap()) == bestRef
                {
                    *(bestHitReads.get_mut(*(ind.get(j).unwrap()) as usize).unwrap()) += 1.0 / numBestRef as f64;

                    if *(xNorm.get(j).unwrap()) >= 0.5
                    {
                        *(level1Reads.get_mut(*(ind.get(j).unwrap()) as usize).unwrap()) += 1.0;
                    }
                    else if *(xNorm.get(j).unwrap()) >= 0.01
                    {
                        *(level2Reads.get_mut(*(ind.get(j).unwrap()) as usize).unwrap()) += 1.0;
                    }
                }
            }

        }

        let readCount = reads.len();

        let bestHit: Vec<f64> = bestHitReads.iter().map(|val| {val.clone() / readCount as f64}).collect();
        let level1: Vec<f64> = level1Reads.iter().map(|val| {val.clone() / readCount as f64}).collect();
        let level2: Vec<f64> = level2Reads.iter().map(|val| {val.clone() / readCount as f64}).collect();

        return (bestHitReads, bestHit, level1, level2);
    }
}

///wrapper for the em function
mod EM
{
    use std::{collections::HashMap};

    pub fn em
    (
        u: &HashMap<i32, (i32, f64)>,
        mut nu: HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        genomes: &Vec<String>,
        maxIter: i32,
        epsilon: f64,
        piPrior: f64,
        thetaPrior: f64
    )
    ->
    (
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>
    )
    {
        let genomeCount = genomes.len();
        let mut pi = vec![1.0 as f64 / genomeCount as f64; genomeCount];
        let mut initPi = pi.clone();
        let mut theta = pi.clone();

        let mut piSum0 = vec![0.0; genomeCount];

        let uWeights: Vec<f64> = u.iter().map(|entry| {(*(entry.1)).1}).collect();
        let mut maxUWeights = 0.0;
        let mut uTotal = 0.0;

        if !uWeights.is_empty()
        {
            maxUWeights = uWeights.iter().cloned().fold(-1./0. /* -inf */, f64::max);
            uTotal = uWeights.iter().sum();
        }

        for i in u.keys()
        {
            piSum0[u.get(i).unwrap().0 as usize] += u.get(i).unwrap().1.clone();
        }
        

        let nuWeights: Vec<f64> = nu.iter().map(|entry| {(*(entry.1)).3}).collect();
        let mut maxNuWeights = 0.0;
        let mut nuTotal = 0.0;

        if !nuWeights.is_empty()
        {
            maxNuWeights = nuWeights.iter().cloned().fold(-1./0. /* -inf */, f64::max);
            nuTotal = nuWeights.iter().sum();
        }

        let priorWeight = f64::max(maxUWeights, maxNuWeights);
        let mut nuLength = nu.len();

        if nuLength == 0
        {
            nuLength = 1;
        }

        //EM iterations
        for i in 0..maxIter
        {
            let piOld = pi.clone();
            let mut thetaSum = vec![0.0; genomeCount];

            //E step
            for j in nu.clone().keys()
            {
                
                let z = nu.get(j).unwrap().clone();

                //A set of any genome mapping with j
                let ind = &z.0;

                //Get relevant pis for the read
                let piTemp: Vec<f64> = ind.iter().map(|val| {pi[*val as usize].clone()}).collect();
                
                //Get relevant thetas for the read
                let thetaTemp: Vec<f64> = ind.iter().map(|val| {theta[*val as usize].clone()}).collect();

                //Calculate non-normalized xs
                let mut xTemp: Vec<f64> = Vec::new();

                for k in 0..ind.len()
                {
                    xTemp.push(piTemp[k] * thetaTemp[k] * z.1[k]);
                }

                let xSum: f64 = xTemp.iter().sum();

                //Avoid dividing by 0 at all times
                let xNorm: Vec<f64>;

                if xSum == 0.0
                {
                    xNorm = vec![0.0; xTemp.len()];
                }
                else
                {
                    xNorm = xTemp.iter().map(|val| {val / xSum}).collect();
                }

                //Update x in nu
                nu.get_mut(j).unwrap().2 = xNorm.clone();

                for (k, _) in ind.iter().enumerate()
                {
                    thetaSum[ind[k] as usize] += xNorm[k] * nu.get(j).unwrap().3;   
                }
            }

            //M step
            let piSum: Vec<f64> = thetaSum.iter().enumerate().map(|(idx, _)| {thetaSum[idx].clone() + piSum0[idx].clone()}).collect();
            let pip = piPrior * priorWeight;

            //update pi
            pi = piSum.iter().map(|val| {((*val as f64) + pip) / (uTotal + nuTotal + (pip * piSum.len() as f64))}).collect();

            if i == 0
            {
                initPi = pi.clone();
            }

            let thetaP = thetaPrior * priorWeight;

            let mut nuTotalDiv = nuTotal.clone();

            if nuTotalDiv == 0 as f64
            {
                nuTotalDiv = 1 as f64;
            }

            theta = thetaSum.iter().map(|val| {(*val + thetaP)/(nuTotalDiv + (thetaP * thetaSum.len() as f64))}).collect();
            

            let mut cutoff = 0.0;

            for (k, _) in pi.iter().enumerate()
            {
                cutoff += (piOld[k] - pi[k]).abs();
            }

            if cutoff <= epsilon || nuLength == 1
            {
                break;
            }
        }

        return (initPi, pi, theta, nu);
    }
}

///Wrapper for the rewriteAlign function and adjacent code
mod RewriteAlign
{
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::LineWriter;
    use std::io::Write;
    use std::io::BufReader;
    use crate::ParseSam::*;

    pub fn rewriteAlign
    (
        u: &HashMap<i32, (i32, f64)>,
        nu: &HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        SAMPath: &str,
        pScoreCutoff: &f64,
        path: &str
    )
    {
        let mut readIdDict: HashMap<String, i32> = HashMap::new();
        let mut refIdDict: HashMap<String, i32> = HashMap::new();

        let mut genomes: Vec<String> = Vec::new();
        let mut read: Vec<String> = Vec::new();

        let mut refCount = 0;
        let mut readCount = 0;

        let oldFile = File::open(SAMPath).expect("Invalid file");
        let mut SAMReader = BufReader::new(oldFile);
        let newFile = File::create(path).expect("unable to create file");
        let mut SAMWriter = LineWriter::new(newFile);

        //for line in parseSam
        loop
        {
            let samLine = parseSAM(&mut SAMReader, Some(*pScoreCutoff));

            let samLine = match samLine
            {
                parseResult::Ok(line) => line,
                parseResult::Ignore => continue,
                parseResult::EOF => break,
                parseResult::Err(_) => panic!("unable to read oldFile in RewriteAlign::rewriteAlign")
            };

            let mut refIndex = refIdDict.get(&samLine.ref_id).unwrap_or(&-1).clone();

            if refIndex == -1
            {
                refIndex = refCount.clone();
                refIdDict.insert(samLine.ref_id.clone(), refIndex);
                genomes.push(samLine.ref_id);
                refCount += 1;
            }

            let mut readIndex = readIdDict.get(&samLine.read_id).unwrap_or(&-1).clone();
            
            if readIndex == -1
            {
                // hold on to this new read
                // first, wrap previous read profile and see if any previous read has a
                // same profile with that!

                readIndex = readCount.clone();
                readIdDict.insert(samLine.read_id.clone(), readIndex);
                read.push(samLine.read_id);
                readCount += 1;

                if u.contains_key(&readIndex)
                {
                    SAMWriter.write(samLine.line.as_bytes()).expect("unable to write to newFile in RewriteAlign::rewriteAlign");
                    continue;
                }
            }

            if nu.contains_key(&readIndex)
            {
                if findUpdatedScore(&nu, readIndex, refIndex) < *pScoreCutoff
                {
                    continue;
                }
                SAMWriter.write(samLine.line.as_bytes()).expect("unable to write to newFile in RewriteAlign::rewriteAlign");
            }
        }
    }


    fn findUpdatedScore
    (
        nu: &HashMap<i32, (Vec<i32>, Vec<f64>, Vec<f64>, f64)>,
        readIndex: i32,
        refIndex: i32
    )
    ->
        f64
    {
        let v1 = match nu.get(&readIndex)
        {
            Some(val) => val,
            None => return 0.0
        };

        let mut idx: usize = 0;

        for (i, el) in v1.0.iter().enumerate()
        {
            if *el == refIndex
            {
                idx = i;
                break;
            }
        }

        return nu.get(&readIndex).unwrap().2.get(idx).unwrap().clone();
    }
}


use pyo3::prelude::*;

#[pymodule]
///pyo3 interface
fn virtool_expectation_maximization(_py: Python, m: &PyModule) -> PyResult<()>
{
    m.add_function(wrap_pyfunction!(run, m)?)?;
    return Ok(());
}


#[pyfunction]
///Entry point for the virtool_expectation_maximization python module
pub fn run
(
    _py: Python,
    SAMPath: String,
    reassignedPath: String,
    pScoreCutoff: f64,
)
->
(
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
    Vec<String>
)
{


    let (u, nu, refs, reads) = BuildMatrix::buildMatrix(SAMPath.as_str(), None);
    
    let (bestHitInitialReads, bestHitInitial, level1Initial, level2Initial) = BestHit::computeBestHit(&u, &nu, &refs, &reads);

    let (initPi, pi, _, nu) = EM::em(&u, nu, &refs, 50, 1e-7, 0.0, 0.0);
    
    let (bestHitFinalReads, bestHitFinal, level1Final, level2Final) = BestHit::computeBestHit(&u, &nu, &refs, &reads);
    
    RewriteAlign::rewriteAlign(&u, &nu, SAMPath.as_str(), &pScoreCutoff, &reassignedPath);

 
    return
    (
        bestHitInitialReads,
        bestHitInitial,
        
        level1Initial,
        level2Initial,
        
        bestHitFinalReads,
        bestHitFinal,

        level1Final,
        level2Final,

        initPi,
        pi,
        
        refs,
        reads
    );
}

///tests and whatnot
#[cfg(test)]
mod tests
{
    #![allow(unused)]
    use std::io::BufReader;
    use std::fs::File;
    use crate::BestHit::*;
    use crate::ParseSam::*;
    use crate::BuildMatrix::*;
    use crate::EM::*;
    use crate::RewriteAlign::*;

    // commented out because the pyo3 interface run
    // accepts different arguments than the rust run

    // #[test]
    // fn testRun()
    // {
    //     let (
    //         bestHitInitialReads,
    //         bestHitInitial,
            
    //         level1Initial,
    //         level2Initial,
            
    //         bestHitFinalReads,
    //         bestHitFinal,
    
    //         level1Final,
    //         level2Final,
    
    //         initPi,
    //         pi,
            
    //         refs,
    //         reads
    //     ) = super::run(String::from("TestFiles/test_al.sam"), String::from("TestFiles/rewrite.sam"), 0.01);

    //     let mut records: Vec<(
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         f64,
    //         String,
    //         String)> = Vec::new();

    //     for i in 0..bestHitFinalReads.len()
    //     {
    //         records.push((
    //             bestHitInitialReads[i],
    //             bestHitInitial[i],
            
    //             level1Initial[i],
    //             level2Initial[i],
            
    //             bestHitFinalReads[i],
    //             bestHitFinal[i],
    
    //             level1Final[i],
    //             level2Final[i],
    
    //             initPi[i],
    //             pi[i],
            
    //             refs[i].clone(),
    //             reads[i].clone()
    //         ));
    //     }

    //     println!("break")
        
    // }

    #[test]
    fn testRewriteAlign()
    {
        let (u, nu, refs, reads) = buildMatrix("TestFiles/test_al.sam", None);
        let (initPi, pi, theta, nu) = em(&u, nu, &refs, 5, 1e-7, 0.0, 0.0);
        rewriteAlign(&u, &nu, "TestFiles/test_al.sam", &0.01, "TestFiles/rewrite.sam");
        println!("pause!");
    }

    ///tests the em function
    #[test]
    fn testem()
    {
        let (u, nu, refs, reads) = buildMatrix("TestFiles/test_al.sam", None);
        let (initPi, pi, theta, nu) = em(&u, nu, &refs, 5, 1e-7, 0.0, 0.0);
        println!("pause!");
    }

    #[test]
    fn testBestHit()
    {
        let (u, nu, refs, reads) = buildMatrix("TestFiles/test_al.sam", None);
        let (bestHitReads, bestHit, level1, level2) = computeBestHit(&u, &nu, &refs, &reads);
        println!("pause");
    }

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