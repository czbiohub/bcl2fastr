//! Represents SampleSheet information as struct with lanes
//! that have index to sample mappings with all indices included within distance 1
//! of original index or distance 0 if overlapping indices are present within distance 1

use std::path::PathBuf;
use std::collections::{HashMap, HashSet};

use itertools::Itertools;


/// SampleData maps from lane number to the index map for the lane
pub type SampleData = HashMap<u32, Lane>;


/// Lane has one or two maps from index to a tuple of (sample name, corrected index)
/// If there is only one index, index2 will be empty and shouldn't be used
#[derive(Debug, PartialEq)]
pub struct Lane {
    index: HashMap<Vec<u8>, (String, String)>,
    index2: HashMap<Vec<u8>, (String, String)>
}


impl Lane {
    /// Look up a sample given a vector of indicies (themselves vectors)
    pub fn get_sample(&self, indices: Vec<Vec<u8>>) -> Option<(String, String)> {
        match indices.len() {
            1 => self.get_1index_sample(&indices[0]),
            2 => self.get_2index_sample(&indices[0], &indices[1]),
            _ => panic!("Got {} indices?!", indices.len())
        }
    }

    /// helper function for when there is one index
    fn get_1index_sample(&self, i: &Vec<u8>) -> Option<(String, String)> {
        match self.index.get(i) {
            Some(sample_t) => return Some(sample_t.clone()),
            None => return None
        };
    }

    /// helper function for when there are two indices
    fn get_2index_sample(&self, i: &Vec<u8>, i2: &Vec<u8>) -> Option<(String, String)> {
        let (index1_sample, og_index) = match self.index.get(i) {
            Some(sample_t) => sample_t,
            None => return None
        };

        let (index2_sample, og_index2) = match self.index2.get(i2) {
            Some(sample_t) => sample_t,
            None => return None
        };

        if index1_sample == index2_sample {
            return Some((index1_sample.clone(), format!("{}+{}", og_index, og_index2)))
        } else {
            return None
        }
    }
}


/// given a set of indices, return a new set of indices which include everything
/// within hamming distance 1
fn hamming_set(index_set: &HashSet<Vec<u8>>) -> HashSet<Vec<u8>> {
    let nucleotides = [b'A', b'C', b'G', b'T', b'N'];

    let new_set: HashSet<_> = index_set.iter().cloned().map(
        |index| {
            let mut this_set = HashSet::new();
            for i in 0..index.len() {
                for c in nucleotides.iter().cloned() {
                    let mut new_index = index.clone();
                    new_index[i] = c;
                    this_set.insert(new_index);
                }
            }
            this_set
        }
    ).flatten().collect();

    new_set
}


/// Function to check for overlaps between the sets of sample indices. If there are
/// two indices, then an overlap between one is allowed as long as the second index
/// is sufficient to distinguish them.
fn check_conflict(
    index_sets: &Vec<HashSet<Vec<u8>>>, index2_sets: &Vec<HashSet<Vec<u8>>>
) -> bool {
    let sample_clash: HashSet<_> = index_sets.iter()
        .enumerate()
        .tuple_combinations()
        .filter_map(
            |((s1, hset1), (s2, hset2))| 
            {
                if hset1.intersection(hset2).count() > 0 {
                    Some((std::cmp::min(s1, s2), std::cmp::max(s1, s2)))
                } else {
                    None
                }
            }
        )
        .collect();

    if index2_sets.len() == 0 && sample_clash.len() > 0 {
        return true
    }

    let sample_clash2: HashSet<_> = index2_sets.iter()
        .enumerate()
        .tuple_combinations()
        .filter_map( 
            |((s1, hset1), (s2, hset2))| {
                if hset1.intersection(hset2).count() > 0 {
                    Some((std::cmp::min(s1, s2), std::cmp::max(s1, s2)))
                } else {
                    None
                }
            }
        )
        .collect();

    if sample_clash.intersection(&sample_clash2).count() > 0 {
        return true
    } else {
        return false
    }
}


/// Function to go from a lane worth of sample and index vectors to a Lane struct
///  which will include the necessary error-correction, up to some limit `max_distance`
fn make_sample_maps(
    sample_names: &Vec<String>,
    index_vec: &Vec<Vec<u8>>,
    index2_vec: &Vec<Vec<u8>>,
    max_distance: usize
) -> Lane {
    // index2_vec should either be full or empty, nothing in between
    assert!(
        index2_vec.len() == index_vec.len() || index2_vec.len() == 0,
        "Samplesheet is missing index2 for some samples"
    );

    // start at distance 0: just map samples to indices
    let mut index_hash_sets: Vec<_> = index_vec.iter()
        .map(
            |index| [index.clone()].iter().cloned().collect::<HashSet<_>>()
        )
        .collect();

    let mut index2_hash_sets: Vec<_> = index2_vec.iter()
        .map(
            |index2| [index2.clone()].iter().cloned().collect::<HashSet<_>>()
        )
        .collect();

    if check_conflict(&index_hash_sets, &index2_hash_sets) {
        panic!("Can't demux two different samples using the same indices");
    }

    for _ in 1..=max_distance {
        let new_index_hash_sets: Vec<_> = index_hash_sets.iter()
            .map(hamming_set)
            .collect();

        let new_index2_hash_sets: Vec<_> = index2_hash_sets.iter()
            .map(hamming_set)
            .collect();

        if check_conflict(&new_index_hash_sets, &new_index2_hash_sets) {
            break
        }

        index_hash_sets = new_index_hash_sets;
        index2_hash_sets = new_index2_hash_sets;
    }

    // converted to this from utf-8 originally, so we know this is safe
    let index_vec: Vec<_> = index_vec.iter()
        .cloned()
        .map( |index| unsafe { String::from_utf8_unchecked(index) })
        .collect();

    // filter_map so that this is empty if there is no second index
    let index2_vec: Vec<_> = index2_vec.iter()
        .cloned()
        .map( |index| unsafe { String::from_utf8_unchecked(index) })
        .collect();

    let mut index_map = HashMap::new();
    let mut index2_map = HashMap::new();

    for ((sample_name, og_index), hash_set) in sample_names.iter()
        .zip(index_vec.iter())
        .zip(index_hash_sets) {
        for index in hash_set.iter().cloned() {
            index_map.insert(index, (sample_name.clone(), og_index.clone()));
        }
    }

    for ((sample_name, og_index2), hash_set) in sample_names.iter()
        .zip(index2_vec.iter())
        .zip(index2_hash_sets) {

        for index2 in hash_set.iter().cloned() {
            index2_map.insert(index2, (sample_name.clone(), og_index2.clone()));
        }
    }

    Lane { index: index_map, index2: index2_map }
}


/// loads a sample sheet and converts it into a SampleData struct. Our version
/// automatically determines the mismatch rate that prevents conflicts, up to
/// a specified maximum
pub fn read_samplesheet(
    samplesheet: PathBuf, max_distance: usize
) -> std::io::Result<SampleData> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_path(samplesheet)?;

    let rows: Vec<_> = rdr.records()
        .filter_map(|r| match r { 
            Ok(r) => Some(r),
            Err(_) => None,
        })
        .skip_while(|r| &r[0] != "[Data]")
        .collect();

    assert!(rows.len() > 2, "No samples found in samplesheet");

    // check for required columns before we start processing
    {
        let row_set: HashSet<_> = rows[1].iter().collect();
        if !row_set.contains(&"Sample_Name") {
            panic!("Samplesheet does not have a Sample_Name column")
        }

        if !row_set.contains(&"Index") {
            panic!("Samplesheet does not have an Index column")
        }
    }

    // collect samples per-lane (or in one big lane if there is no lane column)
    let mut lanes = HashMap::new();

    for record in rows[2..].iter()
        .map( |r| 
            rows[1].iter()
                .zip(r.iter())
                .collect::<HashMap<_, _>>()
        ) {
        let lane: u32 = match record.get(&"Lane") {
            Some(lane) => lane.parse().unwrap(),
            None => 0,
        };
        let (sample_name, sample_idx, sample_idx2) = lanes.entry(lane).or_insert_with(
            || (Vec::new(), Vec::new(), Vec::new())
        );

        sample_name.push(record.get(&"Sample_Name").unwrap().to_string());
        sample_idx.push(record.get(&"Index").unwrap().as_bytes().to_vec());
        match record.get(&"Index2") {
            Some(&index2) => sample_idx2.push(index2.as_bytes().to_vec()),
            None => ()
        };
    }

    let sample_data: HashMap<_, _> = lanes.iter()
        .map( |(&i, (sample_name, sample_idx, sample_idx2))|
            (i, make_sample_maps(sample_name, sample_idx, sample_idx2, max_distance))
        )
        .collect();
 
    Ok(sample_data)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn with_conflict_no_index2_w_lanes() {
        let samplesheet = PathBuf::from(
            "test_data/sample_data/w_conflict_no_index2_w_lanes.csv"
        );
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let test_contents = fs::read_to_string(
            "test_data/hamming_distance_1_test.txt"
        ).unwrap();

        let hamming_set: HashSet<_> = test_contents.split_whitespace()
            .map(|i| i.as_bytes().to_vec())
            .collect();

        let expected_lane1_index: HashMap<_, _> = hamming_set.iter()
            .cloned()
            .map( |ix|  (ix, ("sample_1".to_string(), "ACTGCGAA".to_string())))
            .collect();

        let expected_lane1 = Lane {
            index: expected_lane1_index, index2: HashMap::new()
        };

        let test_contents = fs::read_to_string(
            "test_data/hamming_distance_1_test2.txt"
        ).unwrap();

        let hamming_set2: HashSet<_> = test_contents.split_whitespace()
            .map(|i| i.as_bytes().to_vec())
            .collect();

        let expected_lane2_index: HashMap<_, _> = hamming_set2.iter()
            .cloned()
            .map( |ix|  (ix, ("sample_2".to_string(), "ACTCATCC".to_string())))
            .collect();

        let expected_lane2 = Lane {
            index: expected_lane2_index, index2: HashMap::new()
        };

        let mut expected_sampledata = HashMap::new();

        expected_sampledata.insert(1, expected_lane1);
        expected_sampledata.insert(2, expected_lane2);

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn with_conflict_w_index2_w_lanes() {
        let samplesheet = PathBuf::from(
            "test_data/sample_data/w_conflict_w_index2_w_lanes.csv"
        );
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let index = vec![b"CCCCT".to_vec(), b"CCCCC".to_vec()];
        let index2 = vec![b"AAAAT".to_vec(), b"AAAAA".to_vec()];

        let expected_index: HashMap<_, _> = index.iter()
            .zip(samples.iter().cloned())
            .map( |(ix, sample_name)| 
                (ix.clone(), (sample_name, String::from_utf8(ix.clone()).unwrap()))
            )
            .collect();

        let expected_index2: HashMap<_, _> = index2.iter()
            .zip(samples.iter().cloned())
            .map( |(ix, sample_name)| 
                (ix.clone(), (sample_name, String::from_utf8(ix.clone()).unwrap()))
            )
            .collect();

        let lane_map = sampledata.get(&1).unwrap();

        assert_eq!(lane_map.index, expected_index);
        assert_eq!(lane_map.index2, expected_index2);
    }

    #[test]
    fn with_conflict_in_separate_lanes_w_index2() {
        let samplesheet = PathBuf::from(
            "test_data/sample_data/w_conflict_in_separate_lanes_w_index2.csv"
        );
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];

        let index = "CCCCT".to_string();
        let index_set: Vec<_> = vec![
            "CCTCT", "CGCCT", "CCCCT", "CCCCA", "CCACT", "CNCCT", "GCCCT",
            "CCGCT", "ACCCT", "TCCCT", "CTCCT", "CCCGT", "NCCCT", "CCNCT",
            "CCCTT", "CCCCG", "CACCT", "CCCCN", "CCCNT", "CCCCC", "CCCAT",
        ].iter().map(|i| i.as_bytes().to_vec()).collect();

        let index2 = "AAAAA".to_string();
        let index_set2: Vec<_> = vec![
            "AAAGA", "AGAAA", "AAATA", "AAAAG", "AACAA", "TAAAA", "AAANA",
            "CAAAA", "NAAAA", "AANAA", "ATAAA", "AAACA", "ACAAA", "AAAAN",
            "GAAAA", "AATAA", "AAGAA", "AAAAA", "AAAAT", "AAAAC", "ANAAA"
        ].iter().map(|i| i.as_bytes().to_vec()).collect();


        let expected_index: Vec<HashMap<_, _>> = samples.iter()
            .map( |sample| {
                index_set.iter()
                    .cloned()
                    .map( |ix|  (ix, (sample.clone(), index.clone())))
                    .collect()
            }).collect();

        let expected_index2: Vec<HashMap<_, _>> = samples.iter()
            .map( |sample| {
                index_set2.iter()
                    .cloned()
                    .map( |ix|  (ix, (sample.clone(), index2.clone())))
                    .collect()
            }).collect();

        assert_eq!(sampledata.len(), 2);
        for (lane, lane_map) in sampledata {
            assert_eq!(lane_map.index, expected_index[(lane - 1) as usize]);
            assert_eq!(lane_map.index2, expected_index2[(lane - 1) as usize]);
        }
    }

    #[test]
    fn no_conflict_w_index2() {
        let samplesheet = PathBuf::from(
            "test_data/sample_data/no_conflict_w_index2.csv"
        );
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];

        let mut expected_index = HashMap::new();

        for (k, v) in vec![
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC", "NGGGG",
            "GGGCG", "GGNGG", "GCGGG", "GNGGG", "GGGNG", "GGGGA", "GGGGN",
            "GAGGG", "GGGTG", "GGGGG", "GGCGG", "CGGGG", "GTGGG", "GGGGT"
        ].iter().map(
            |i| (i.as_bytes().to_vec(), (sample_names[0].clone(), "GGGGG".to_string()))
        ) {
            expected_index.insert(k, v);
        }

        for (k, v) in vec![
            "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT", "TTTTT",
            "NTTTT", "ATTTT", "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT",
            "TTTTA", "TTTTC", "TTTTG", "CTTTT", "TTTGT", "TTTTN", "TTATT",
        ].iter().map(
            |i| (i.as_bytes().to_vec(), (sample_names[1].clone(), "TTTTT".to_string()))
        ) {
            expected_index.insert(k, v);
        }

        let mut expected_index2 = HashMap::new();

        for (k, v) in vec![
            "ATAAA", "CAAAA", "AAGAA", "AACAA", "AAAAC", "AAAAG", "AAAAA",
            "AAATA", "ANAAA", "AAAGA", "AGAAA", "GAAAA", "NAAAA", "ACAAA",
            "AANAA", "AAAAN", "AAAAT", "AAACA", "AAANA", "TAAAA", "AATAA"   
        ].iter().map(
            |i| (i.as_bytes().to_vec(), (sample_names[0].clone(), "AAAAA".to_string()))
        ) {
            expected_index2.insert(k, v);
        }

        for (k, v) in vec![
            "CCGCC", "CCNCC", "CCTCC", "TCCCC", "CGCCC", "CCCNC", "ACCCC",
            "CCCTC", "CCCGC", "CCCAC", "CCCCA", "CACCC", "GCCCC", "CCCCC",
            "CCCCT", "CTCCC", "CNCCC", "NCCCC", "CCCCG", "CCACC", "CCCCN"            
        ].iter().map(
            |i| (i.as_bytes().to_vec(), (sample_names[1].clone(), "CCCCC".to_string()))
        ) {
            expected_index2.insert(k, v);
        }

        let mut expected_sampledata = SampleData::new();

        expected_sampledata.insert(
            0,
            Lane { index: expected_index, index2: expected_index2 }
        );

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn sample_lookup() {
        let samplesheet = PathBuf::from(
            "test_data/sample_data/no_conflict_w_index2.csv"
        );
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();
        let lane = sampledata.get(&0).unwrap();

        assert_eq!(
            lane.get_sample(vec![
                vec![71, 84, 71, 71, 71], vec![65, 65, 65, 65, 65] 
            ]).unwrap(), 
            ("sample_1".to_string(), "GGGGG+AAAAA".to_string())
        );

        assert_eq!(
            lane.get_sample(vec![
                vec![65, 65, 65, 65, 65], vec![71, 84, 71, 71, 71],
            ]),
            None
        );

        assert_eq!(
            lane.get_sample(vec![vec![71, 84, 71, 71, 71]]).unwrap(), 
            ("sample_1".to_string(), "GGGGG".to_string())
        );

        assert_eq!(
            lane.get_sample(vec![vec![65, 65, 65, 65, 65]]),
            None
        );
    }

    #[test]
    fn make_sample_maps() {
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let mut expected_index = HashMap::new();

        for (k, v) in vec![
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC", "NGGGG",
            "GGGCG", "GGNGG", "GCGGG", "GNGGG", "GGGNG", "GGGGA", "GGGGN",
            "GAGGG", "GGGTG", "GGGGG", "GGCGG", "CGGGG", "GTGGG", "GGGGT"
        ].iter().map(
            |i| (i.as_bytes().to_vec(), (sample_names[0].clone(), "GGGGG".to_string()))
        ) {
            expected_index.insert(k, v);
        }

        for (k, v) in vec![
            "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT", "TTTTT",
            "NTTTT", "ATTTT", "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT",
            "TTTTA", "TTTTC", "TTTTG", "CTTTT", "TTTGT", "TTTTN", "TTATT"
        ].iter().map(
            |i| (i.as_bytes().to_vec(), (sample_names[1].clone(), "TTTTT".to_string()))
        ) {
            expected_index.insert(k, v);
        }

        let index_vec = vec![b"GGGGG".to_vec(), b"TTTTT".to_vec()];

        let actual_mapping = super::make_sample_maps(
            &sample_names, &index_vec, &vec![], 1
        );

        assert_eq!(actual_mapping.index, expected_index);
    }

    #[test]
    fn make_sample_maps_conflict() {
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let index_vec = vec![b"ACTG".to_vec(), b"ACTC".to_vec()];

        let expected_index: HashMap<_, _> = sample_names.iter()
            .zip(index_vec.iter().cloned())
            .map( |(sample, ix)| 
                (ix.clone(), (sample.clone(), String::from_utf8(ix.clone()).unwrap()))
            )
            .collect();

        let actual_mapping = super::make_sample_maps(
            &sample_names, &index_vec, &vec![], 1
        );

        assert_eq!(actual_mapping.index, expected_index);
    }

    #[test]
    fn hamming_set_distance_1() {
        let initial_set: HashSet<_> = [b"ACTGCGAA".to_vec()].iter().cloned().collect();

        let actual_hammingset = hamming_set(&initial_set);

        let test_contents = fs::read_to_string(
            "test_data/hamming_distance_1_test.txt"
        ).unwrap();

        let expected_hammingset: HashSet<_> = test_contents.split_whitespace()
            .map(|i| i.as_bytes().to_vec())
            .collect();

        assert_eq!(actual_hammingset, expected_hammingset);
    }

    #[test]
    fn hamming_set_distance_2() {
        let initial_set: HashSet<_> = [b"ACTGCGAA".to_vec()].iter().cloned().collect();

        let actual_hammingset = hamming_set(&hamming_set(&initial_set));

        let test_contents = fs::read_to_string(
            "test_data/hamming_distance_2_test.txt"
        ).unwrap();

        let expected_hammingset: HashSet<_> = test_contents.split_whitespace()
            .map(|i| i.as_bytes().to_vec())
            .collect();

        assert_eq!(actual_hammingset, expected_hammingset);
    }

    #[test]
    #[should_panic(
        expected = r#"No such file or directory"#
    )]
    fn no_file() {
        let samplesheet = PathBuf::from("test_data/sample_data/no_file.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"No samples found in samplesheet"#
    )]
    fn bad_file() {
        let samplesheet = PathBuf::from("test_data/sample_data/bad_file.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"No samples found in samplesheet"#
    )]
    fn empty_file() {
        let samplesheet = PathBuf::from("test_data/empty_file");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Can't demux two different samples using the same indices"#
    )]
    fn sample_collision() {
        let samplesheet = PathBuf::from("test_data/sample_data/sample_collision.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Samplesheet does not have a Sample_Name column"#
    )]
    fn no_sample() {
        let samplesheet = PathBuf::from("test_data/sample_data/no_sample_name.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Samplesheet does not have an Index column"#
    )]
    fn no_index() {
        let samplesheet = PathBuf::from("test_data/sample_data/no_index.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }
}
