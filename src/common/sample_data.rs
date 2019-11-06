//! Represents SampleSheet information as struct with lanes
//! that have index to sample mappings with all indices included within distance 1
//! of original index or distance 0 if overlapping indices are present within distance 1


use std::vec::Vec;
use std::path::PathBuf;
use std::collections::HashMap;
use std::collections::HashSet;

use itertools::{Itertools, iproduct};
use csv::Reader;


fn product_repeat(seq: &Vec<char>, k: usize) -> Vec<String> {
    // return all possible substitution number of variations with number of replacements k
    // :param seq: vector of unique characters in alphabet
    // :param k: integer value of substitutions to be made
    // :return: Vector of characters with all possible combinations given k
    match k {
        0 => vec![],
        1 => (0..&seq.len() -1).map(|c| c.to_string()).collect(),
        2 => iproduct!((0..&seq.len() -1), (0..&seq.len() -1)).map(|(a, b)| format!("{}{}", a, b)).collect(),
        _ => iproduct!(product_repeat(&seq.clone(), k - 1), (0..&seq.len() -1)).map(|(a, b)| format!("{}{}", a, b)).collect(),
    }
}


fn hamming_set(index: &str, distance: usize) -> HashSet<String> {
    // Given an index of bases in {ACGTN}, generate all indices within hamming
    // distance of the input
    // :param index: string representing the index sequence
    // :param distance: maximum distance to allow
    // :param include_N: include N when generating possible indexes
    // :return: HashSet of indexes within hamming distance
   
    let nucleotides: Vec<char> = String::from("ACTGN")
        .chars()
        .collect();
    
    let mut hamming_set = HashSet::new();
    // add original index to include hamming distance of 0
    hamming_set
        .insert(
            index
                .to_string()
        );

    for dist in 0..distance+1{
        // create combinations of all indexes to change in original string
        let combinations  = (0..index.len()).combinations(dist);
        for positions in combinations {
            // get all replacements available with given alphabet
            let cartesian_product = product_repeat(&nucleotides, dist);
            for p in cartesian_product {
                let replacements: Vec<i32> = p
                    .chars()
                    .map(
                        |c| c
                            .to_string()
                            .parse()
                            .unwrap()
                    ).collect();
                
                let index_vec: Vec<char> = index
                    .chars()
                    .collect();
                
                let mut index_hamming = HashMap::new();
                
                // create hashmap of original index to update chars later on
                for (pos, e) in index_vec
                    .iter()
                    .enumerate() {
                        index_hamming.insert(pos, e);
                }

                for (p, r) in (
                    positions
                        .iter()
                ).zip(replacements) {
                    let current_val = **index_hamming.get(&p).unwrap();
                    if current_val == nucleotides[*&r as usize] {
                        index_hamming
                            .insert(
                                *p,
                                nucleotides
                                    .last()
                                    .unwrap()
                            );
                    } else {
                        index_hamming.insert(*p, &nucleotides[*&r as usize]);
                    }
                }
                // re-assemble new string with hamming distance from hashmap
                let mut index_hamming_str = String::new();
                // order keys and create distance index
                for a in 0..index_hamming.keys().len() {
                    index_hamming_str.push(*index_hamming[&a])
                }
                hamming_set.insert(index_hamming_str);
            }
        }
    }
    hamming_set
}


fn index_sample_map(indices: &Vec<String>, samples: &Vec<String>, max_dist: usize) -> Result<HashMap<String, String>, HashMap<String, String>> {

    // Given a sequence of indices and samples, return mapping of all possible indices within distance max_dist
    // of the original index to each corresponding sample. if overlapping indices exist, return index to sample mapping
    // with max_dist of 0.
    // :param indices: vector of indices (within a lane)
    // :param samples: vector of samples (within a lane)    
    // :param max_dist: maximum distance to allow
    // :return:
    //     Ok => HashMap of all indices within max_dist to sample
    //     Err => HashMap of all indices within distance 0 to sample
    
    let mut all_hsets = HashMap::new();
    for d in 0..max_dist+1 {

        let h_sets: Vec<HashSet<String>> = indices
            .iter()
            .map(
                |x| hamming_set(x, d)
            ).collect();
        
        // samples to indices within d
        let sample_to_hsets: Vec<(String, HashSet<String>)> = samples
            .iter()
            .cloned()
            .zip(
                h_sets.iter()
                    .cloned()
            ).collect();
        
        all_hsets.entry(d).or_insert(sample_to_hsets.clone());
        
        for (hset1, hset2) in h_sets
            .iter()
            .flatten()
            .tuple_combinations() {            
                // check for all pairs if any equal
                if hset1 == hset2 { 
                    let mapping: HashMap<String, String> = indices
                        .iter()
                        .cloned()
                        .zip(
                            samples.iter()
                                .cloned()
                        ).collect();
                    return Err(mapping);
            }  
        }
    }
    
    let mut mapping = HashMap::new();
    let highest_distance = all_hsets.get(&max_dist).unwrap();
    for (sample, indices) in highest_distance {
        for i in indices.iter().cloned() {
            mapping.insert(i, sample.clone());
        }
    }
    Ok(mapping)
}


// Represents samplesheet data with lanes
#[derive(Debug)]
pub struct SampleData {
    lanes: Vec<Lane> 
}

// Represents Lane where each lane has a specific index to sample mapping
#[derive(Debug)]
pub struct Lane  {
    // lane id
    lane: u32,
    // how far of a distance can an index be away from the orignal index
    // to be considered belonging to a sample
    hamming_distance: usize,
    // mapping of index to sample name for demuxing
    index_to_sample: HashMap<String, String>
}


impl SampleData {
// loads samplesheet for NovaSeqRun given path
    pub fn read_path(samplesheet: PathBuf) -> std::io::Result<SampleData> {
        let mut rdr = Reader::from_path(samplesheet)?;
        let mut rows = rdr.records();
        let row_header = rows
            .next()
            .unwrap()?;
        
        let mut row_map = HashMap::new();
        for (pos, e) in row_header
            .iter()
            .enumerate() {
	        row_map.insert(
                    e, pos
                );
        }

        let mut lanes = HashMap::new();

        for row in rows {
            let record = row?;

            // check if split lanes
            let lane: u32 = if row_map.contains_key(&"Lane") {
                record[*row_map.get(&"Lane").unwrap()]
                    .parse()
                    .unwrap()               
            } else {
            // no split lanes
                1
            };
            
            let lane_entry = lanes
                .entry(lane)
                .or_insert(
                    HashMap::new()
                );
            
            let sample_name = String::from(
                &record[*row_map.get(&"Sample_Name").unwrap()]
            );
                
            let sample_idx = String::from(
                &record[*row_map.get(&"Index").unwrap()]
            );
            
            // push indices and sample names to vec to map later
            let index_entry = lane_entry
                .entry(
                    "indices".to_owned()
                ).or_insert(Vec::new());
            index_entry.push(sample_idx.clone());

            let samples_entry = lane_entry
                .entry(
                    "samples".to_owned()
                ).or_insert(Vec::new());
            samples_entry.push(sample_name.to_string());

            // check if index2
            if row_map.contains_key(&"Index2") {
                let sample_idx2 = String::from(
                    &record[*row_map.get(&"Index2").unwrap()]
                );
                
                lanes.get_mut(&lane)
                    .unwrap()
                    .get_mut("indices")
                    .unwrap()
                    .push(sample_idx2.clone());
                
                // push sample name again if index2
                lanes.get_mut(&lane)
                    .unwrap()
                    .get_mut("samples")
                    .unwrap()
                    .push(sample_name.to_string());                
            }
        }

        let mut total_lanes = Vec::new();
        for (lane, entries) in lanes {
            let samples = entries.get("samples").unwrap();
            let indices = entries.get("indices").unwrap();

            // check that all indices are far away enough (default distance of 1)
            let hc = index_sample_map(&indices, &samples, 1); 
            
            let (index_to_sample, hamming_distance) = match hc {
                Ok(v) => (v, 1),
                // if conflict then return 1:1 mapping and distance 0
                Err(e) => (e, 0)
            };
            let lane = Lane {
                lane: lane,
                hamming_distance: hamming_distance,
                index_to_sample: index_to_sample
            };
            total_lanes.push(lane);
        }
        let sample = SampleData {
            lanes: total_lanes
        };
        Ok(sample)
    }

}



#[cfg(test)]
mod tests {
    use std::fs;
    use super::*;

    #[test]
    fn test_sampledata_w_conflict_no_index2_w_splitlanes() {
        let samplesheet = PathBuf::from("test_data/sample_data/test_samplesheet_w_conflict_no_index2_w_splitlanes.csv");
        let sampledata = SampleData::read_path(samplesheet).unwrap();
        assert_eq!(sampledata.lanes.len(), 2);
        
    }

    #[test]
    fn test_sampledata_w_conflict_w_index2_w_splitlanes() {
        let samplesheet = PathBuf::from("test_data/sample_data/test_samplesheet_w_conflict_w_index2_w_splitlanes.csv");
        let sampledata = SampleData::read_path(samplesheet).unwrap();
        
        let samples: Vec<String> = vec!("800282", "800282", "800283", "800283").iter().map(|i| i.to_string()).collect();
        let indices: Vec<String> = vec!("CCCCT","TTTTT", "CCCCC", "AAAAA").iter().map(|i| i.to_string()).collect();

        let mut expected_mapping = HashMap::new();
        for (s, i) in samples.iter().zip(indices.iter()) {
            expected_mapping.insert(i.to_owned(), s.to_owned());
        }
        for lane in sampledata.lanes {
            assert_eq!(lane.hamming_distance, 0);
                if lane.lane == 1 {
                   assert_eq!(lane.index_to_sample, expected_mapping); 
                }
        }
    }

    #[test]
    fn test_sampledata_w_conflict_in_seperate_lanes_w_index2() {
        let samplesheet = PathBuf::from("test_data/sample_data/test_samplesheet_w_conflict_in_seperate_lanes_w_index2.csv");
        let sampledata = SampleData::read_path(samplesheet).unwrap();

        let index_set1: HashSet<std::string::String> = vec!(
            "CCTCT", "CGCCT", "CCCCT", "CCCCA", "CCACT", "CNCCT", "GCCCT",
            "CCGCT", "ACCCT", "TCCCT", "CTCCT", "CCCGT", "NCCCT", "CCNCT",
            "CCCTT", "CCCCG", "CACCT", "CCCCN", "CCCNT", "CCCCC", "CCCAT",
            "TATTT", "TTTNT", "TTTTC", "TTTTA", "NTTTT", "TTGTT", "TTTAT",
            "TTCTT", "GTTTT", "TGTTT", "CTTTT", "TTNTT", "TTTTN", "TTTTT",
            "TTATT", "ATTTT", "TTTGT", "TTTTG", "TTTCT", "TCTTT", "TNTTT"
        ).iter().map(|i| i.to_string()).collect();
        
        let index_set2: HashSet<std::string::String> = vec!(
            "CCCCT", "CCACC", "CCCNC", "CCTCC", "CCCGC", "CCCCG", "CTCCC",
            "TCCCC", "ACCCC", "CCCCN", "CCCTC", "NCCCC", "CCCCC", "CCCCA",
            "GCCCC", "CCNCC", "CNCCC", "CGCCC", "CCCAC", "CACCC", "CCGCC",
            "AAAGA", "AGAAA", "AAATA", "AAAAG", "AACAA", "TAAAA", "AAANA",
            "CAAAA", "NAAAA", "AANAA", "ATAAA", "AAACA", "ACAAA", "AAAAN",
            "GAAAA", "AATAA", "AAGAA", "AAAAA", "AAAAT", "AAAAC", "ANAAA"
        ).iter().map(|i| i.to_string()).collect();

        let mut expected_mapping_lane1 = HashMap::new();
        let all_s1: Vec<String> = vec!("800282".to_string(); index_set1.len());
        for (i, ss) in index_set1.iter().zip(all_s1.iter()) {
            expected_mapping_lane1.insert(i.to_owned(), ss.to_owned());
        }

        let mut expected_mapping_lane2 = HashMap::new();
        let all_s2: Vec<String> = vec!("800283".to_string(); index_set2.len());
        for (i, ss) in index_set2.iter().zip(all_s2.iter()) {
            expected_mapping_lane2.insert(i.to_owned(), ss.to_owned());
        }

        assert_eq!(sampledata.lanes.len(), 2);
        for lane in sampledata.lanes {
            assert_eq!(lane.hamming_distance, 1);
            if lane.lane == 1 {
                assert_eq!(lane.index_to_sample, expected_mapping_lane1);
            }
            if lane.lane == 2 {
                assert_eq!(lane.index_to_sample, expected_mapping_lane2);
            }
        }
    }      

    #[test]
    fn test_sampledata_no_conflict_w_index2_no_splitlanes() {
        let samplesheet = PathBuf::from("test_data/sample_data/test_samplesheet_no_conflict_w_index2_no_splitlanes.csv");
        let sampledata = SampleData::read_path(samplesheet).unwrap();

        let samples: Vec<String> = vec!("800282", "800283").iter().map(|i| i.to_string()).collect();
        let index_set1: HashSet<std::string::String> = vec!(
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC","NGGGG",
            "GGGCG", "GGNGG", "GCGGG", "GNGGG", "GGGNG","GGGGA", "GGGGN",
            "GAGGG", "GGGTG", "GGGGG", "GGCGG", "CGGGG", "GTGGG", "GGGGT",
            "ATAAA", "CAAAA", "AAGAA", "AACAA", "AAAAC", "AAAAG", "AAAAA",
            "AAATA", "ANAAA", "AAAGA", "AGAAA", "GAAAA", "NAAAA", "ACAAA",
            "AANAA", "AAAAN", "AAAAT", "AAACA", "AAANA", "TAAAA", "AATAA"
            
        ).iter().map(|i| i.to_string()).collect();
        let index_set2: HashSet<std::string::String> = vec!(
            "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT", "TTTTT",
            "NTTTT", "ATTTT", "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT",
            "TTTTA", "TTTTC", "TTTTG", "CTTTT", "TTTGT", "TTTTN", "TTATT",
            "CCGCC", "CCNCC", "CCTCC", "TCCCC", "CGCCC", "CCCNC", "ACCCC",
            "CCCTC", "CCCGC", "CCCAC", "CCCCA", "CACCC", "GCCCC", "CCCCC",
            "CCCCT", "CTCCC", "CNCCC", "NCCCC", "CCCCG", "CCACC", "CCCCN"            
        ).iter().map(|i| i.to_string()).collect();

        let index_sets = vec!(index_set1, index_set2);
        let mut expected_mapping = HashMap::new();
        for (s, hs) in samples.iter().zip(index_sets.iter()) {
            let all_s: Vec<String> = vec!(s.to_string(); hs.len());
            for (i, ss) in hs.iter().zip(all_s.iter()) {
                expected_mapping.insert(i.to_owned(), ss.to_owned());
            }
        }

        assert_eq!(sampledata.lanes.len(), 1);
        for lane in sampledata.lanes {
            assert_eq!(lane.hamming_distance, 1);
            if lane.lane == 1 {
                assert_eq!(lane.index_to_sample, expected_mapping);
            }
        }
    }      
    
    #[test]
    fn testindex_sample_map_ok() {
        let samples: Vec<String> = vec!("800282", "800283").iter().map(|i| i.to_string()).collect();
        let indices: Vec<String> = vec!("GGGGG", "TTTTT").iter().map(|i| i.to_string()).collect();
        let index_set1: HashSet<std::string::String> = vec!(
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC","NGGGG",
            "GGGCG", "GGNGG", "GCGGG", "GNGGG", "GGGNG","GGGGA", "GGGGN",
            "GAGGG", "GGGTG", "GGGGG", "GGCGG", "CGGGG", "GTGGG", "GGGGT"
        ).iter().map(|i| i.to_string()).collect();
        let index_set2: HashSet<std::string::String> = vec!(
            "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT", "TTTTT",
            "NTTTT", "ATTTT", "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT",
            "TTTTA", "TTTTC", "TTTTG", "CTTTT", "TTTGT", "TTTTN", "TTATT"
        ).iter().map(|i| i.to_string()).collect();

        let index_sets = vec!(index_set1, index_set2);
        let actual_mapping = index_sample_map(&indices, &samples, 1);
        let mut expected_mapping = HashMap::new();
        for (s, hs) in samples.iter().zip(index_sets.iter()) {
            let all_s: Vec<String> = vec!(s.to_string(); hs.len());
            for (i, ss) in hs.iter().zip(all_s.iter()) {
                expected_mapping.insert(i.to_owned(), ss.to_owned());
            }
        }
        
        assert_eq!(actual_mapping.unwrap(), expected_mapping);        
    }

    #[test]
    fn testindex_sample_map_error() {
        let samples: Vec<String> = vec!("800282", "800283").iter().map(|i| i.to_string()).collect();
        let indices: Vec<String> = vec!("ACTG", "ACTC").iter().map(|i| i.to_string()).collect();
        let actual_mapping = index_sample_map(&indices, &samples, 2);
        let expected_mapping: HashMap<String, String> = indices.iter().cloned().zip(samples.iter().cloned()).collect(); 
        assert_eq!(actual_mapping, Err(expected_mapping));   
        
    }
   
    #[test]
    fn testhamming_set_distance_1() {
        // hamming set distance 1
        let actual_hammingset = hamming_set(&String::from("ACTGCGAA"), 1);
        let test_contents = fs::read_to_string("test_data/hamming_distance_1_test.txt").unwrap();
        let expected_hammingset: HashSet<std::string::String> = test_contents.split_whitespace().map(|i| i.to_string()).collect();
        assert_eq!(actual_hammingset, expected_hammingset);        
    }    
    
    #[test]
    fn testhamming_set_distance_2() {
        // hamming set distance 1
        let actual_hammingset = hamming_set(&String::from("ACTGCGAA"), 2);
        let test_contents = fs::read_to_string("test_data/hamming_distance_2_test.txt").unwrap();
        let expected_hammingset: HashSet<std::string::String> = test_contents.split_whitespace().map(|i| i.to_string()).collect();
        assert_eq!(actual_hammingset, expected_hammingset);
    }
        
           
}
