//! Represents SampleSheet information as struct with lanes
//! that have index to sample mappings with all indices included within distance 1
//! of original index or distance 0 if overlapping indices are present within distance 1

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use log::warn;
use ndarray::ArrayView1;
use rayon::prelude::*;

use crate::hamming_set::{check_conflict, hamming_set};

/// SampleData maps from lane number to the index maps for the lane. The maps are
/// chunked into different pieces, each corresponding to a set of samples that will be
/// processed together
pub type SampleData = HashMap<usize, Samples>;

/// The Samples struct has one or two maps that go from potential indices to sample indices.
/// To save space, we save the original names as a vector and use integers everywhere else.
///
/// If there is only one index, the index2 values will be empty.
#[derive(Debug, PartialEq)]
pub struct Samples {
    /// vector of original sample names
    pub sample_names: Vec<String>,
    /// project name for each sample, or None if project wasn't in sample sheet
    pub project_names: Vec<Option<String>>,
    /// vector of exact indices for each sample, for reporting exact match numbers. Contains
    /// a HashSet for each sample because a single sample can have multiple indices listed
    index_vec: Vec<HashSet<Vec<u8>>>,
    /// map from all valid indices (up to error tolerance) to a usize which indexes into the
    /// sample_names vector
    index_map: HashMap<Vec<u8>, HashSet<usize>>,
    /// the same vector of index sets, but for the second index (if used)
    index2_vec: Vec<HashSet<Vec<u8>>>,
    /// another map from valid index to sample usize, for the second index
    index2_map: HashMap<Vec<u8>, HashSet<usize>>,
}

impl Samples {
    /// Look up a sample given a vector of indices
    pub fn get_sample(&self, indices: &[ArrayView1<u8>]) -> Option<(usize, bool)> {
        match indices.len() {
            1 => self.get_1index_sample(indices[0]),
            2 => self.get_2index_sample(indices[0], indices[1]),
            x => panic!("Got {} indices?!", x),
        }
    }

    /// Checks if the indices match any of the samples
    pub fn is_any_sample(&self, indices: &[Vec<u8>]) -> bool {
        match indices.len() {
            1 => self.index_map.contains_key(&indices[0]),
            2 => {
                match (
                    self.index_map.get(&indices[0]),
                    self.index2_map.get(&indices[1]),
                ) {
                    (Some(ix_set), Some(ix2_set)) if !ix_set.is_disjoint(&ix2_set) => true,
                    _ => false,
                }
            }
            x => panic!("Got {} indices?!", x),
        }
    }

    /// helper function for when there is one index. By construction every set will contain
    /// only a single usize element
    fn get_1index_sample(&self, idx: ArrayView1<u8>) -> Option<(usize, bool)> {
        if let Some(ix_set) = self.index_map.get(idx.as_slice().unwrap()) {
            let ix = ix_set.iter().nth(0).unwrap().clone();
            return Some((ix, self.index_vec[ix].contains(idx.as_slice().unwrap())));
        };
        return None;
    }

    /// helper function for when there are two indices. By construction any intersection will
    /// either be empty or contain one usize element. This has a potential "bug" in that it will
    /// say there's an exact match when the two indices match the same sample exactly but were not
    /// specified as a pair in the sample sheet
    fn get_2index_sample(
        &self,
        idx: ArrayView1<u8>,
        idx2: ArrayView1<u8>,
    ) -> Option<(usize, bool)> {
        match (
            self.index_map.get(idx.as_slice().unwrap()),
            self.index2_map.get(idx2.as_slice().unwrap()),
        ) {
            (Some(ix_set), Some(ix2_set)) if !ix_set.is_disjoint(ix2_set) => {
                let ix = ix_set.intersection(ix2_set).nth(0).unwrap().clone();
                let is_exact = self.index_vec[ix].contains(idx.as_slice().unwrap())
                    && self.index2_vec[ix].contains(idx2.as_slice().unwrap());

                Some((ix, is_exact))
            }
            _ => None,
        }
    }
}

/// Function to go from a lane worth of sample and index vectors to a Lane struct
/// which will include the necessary error-correction, up to some limit `max_distance`
fn make_sample_maps(
    sample_names: &[String],
    project_names: &[Option<String>],
    index_vec: &[Vec<u8>],
    index2_vec: &[Vec<u8>],
    max_distance: usize,
) -> Samples {
    // index_vec should be full
    assert_eq!(
        sample_names.len(),
        index_vec.len(),
        "Missing indexes for some samples"
    );

    // index2_vec should either be full or empty, nothing in between
    assert!(
        index2_vec.len() == index_vec.len() || index2_vec.len() == 0,
        "Samplesheet is missing index2 for some samples"
    );

    // same for sample_project: full or empty, nothing in between
    let n_project_names = project_names.iter().cloned().filter_map(|n| n).count();
    assert!(
        n_project_names == sample_names.len() || n_project_names == 0,
        "Samplesheet is missing project names for some samples"
    );

    let mut unique_samples: Vec<_> = sample_names
        .iter()
        .cloned()
        .zip(project_names.iter().cloned())
        .collect();
    unique_samples.sort();
    unique_samples.dedup();

    let (unique_sample_names, unique_project_names): (Vec<_>, Vec<_>) =
        unique_samples.drain(..).unzip();

    let mut index_hash_sets: HashMap<_, _> = HashMap::new();
    let mut index2_hash_sets: HashMap<_, _> = HashMap::new();

    // for each sample_name we combine all the given indexes. Allows for multiple entries
    // to demux into one file
    for (sample_name, ix) in sample_names.iter().zip(index_vec) {
        index_hash_sets
            .entry(sample_name.clone())
            .or_insert(HashSet::new())
            .insert(ix.clone());
    }

    for (sample_name, ix2) in sample_names.iter().zip(index2_vec) {
        index2_hash_sets
            .entry(sample_name.clone())
            .or_insert(HashSet::new())
            .insert(ix2.clone());
    }

    // first check for conflicts at distance 0 (i.e. using same index twice)
    if check_conflict(&index_hash_sets, &index2_hash_sets) {
        panic!("Can't demux two different samples using the same indices");
    }

    let index_vec: Vec<_> = unique_sample_names
        .iter()
        .map(|s| index_hash_sets.get(s).unwrap().clone())
        .collect();
    let index2_vec: Vec<_> = if index2_hash_sets.len() > 0 {
        unique_sample_names
            .iter()
            .map(|s| index2_hash_sets.get(s).unwrap().clone())
            .collect()
    } else {
        Vec::new()
    };

    for i in 1..=max_distance {
        let new_index_hash_sets: HashMap<_, _> = index_hash_sets
            .par_iter()
            .map(|(s, ixs)| (s.clone(), hamming_set(ixs)))
            .collect();
        let new_index2_hash_sets: HashMap<_, _> = index2_hash_sets
            .par_iter()
            .map(|(s, ix2s)| (s.clone(), hamming_set(ix2s)))
            .collect();

        if check_conflict(&new_index_hash_sets, &new_index2_hash_sets) {
            warn!(
                "Warning: conflict at distance {}, using {} instead",
                i,
                i - 1
            );
            break;
        }

        index_hash_sets = new_index_hash_sets;
        index2_hash_sets = new_index2_hash_sets;
    }

    let mut index_map = HashMap::new();
    let mut index2_map = HashMap::new();

    for (i, sample_name) in unique_sample_names.iter().cloned().enumerate() {
        for ix in index_hash_sets.get(&sample_name).unwrap().iter().cloned() {
            index_map.entry(ix).or_insert(HashSet::new()).insert(i);
        }

        if index2_hash_sets.len() > 0 {
            for ix2 in index2_hash_sets.get(&sample_name).unwrap().iter().cloned() {
                index2_map.entry(ix2).or_insert(HashSet::new()).insert(i);
            }
        }
    }

    Samples {
        sample_names: unique_sample_names,
        project_names: unique_project_names,
        index_vec: index_vec,
        index_map: index_map,
        index2_vec: index2_vec,
        index2_map: index2_map,
    }
}

/// loads a sample sheet and converts it into a SampleData struct. Our version
/// automatically determines the mismatch rate that prevents conflicts, up to
/// a specified maximum
pub fn read_samplesheet(samplesheet: PathBuf, max_distance: usize) -> std::io::Result<SampleData> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_path(samplesheet)?;

    // ignore any rows before the Data section
    let rows: Vec<_> = rdr
        .records()
        .filter_map(|r| match r {
            Ok(r) => Some(r),
            Err(e) => panic!("{}", e),
        })
        .skip_while(|r| &r[0] != "[Data]")
        .collect();

    assert!(rows.len() > 2, "No samples found in samplesheet");

    // check for required columns before we start processing
    {
        let row_set: Vec<_> = rows[1].iter().collect();
        if !row_set.contains(&"Sample_Name") {
            panic!("Samplesheet does not have a Sample_Name column")
        }

        if !row_set.contains(&"Index") {
            panic!("Samplesheet does not have an Index column")
        }
    }

    // collect samples per-lane (or in one big lane if there is no lane column)
    let mut lanes = HashMap::new();

    for record in rows[2..]
        .iter()
        .map(|r| rows[1].iter().zip(r.iter()).collect::<HashMap<_, _>>())
    {
        let lane: usize = match record.get(&"Lane") {
            Some(lane) => lane.parse().unwrap(),
            None => 0,
        };

        let (sample_names, project_names, sample_idx, sample_idx2) = lanes
            .entry(lane)
            .or_insert_with(|| (Vec::new(), Vec::new(), Vec::new(), Vec::new()));

        sample_names.push(record.get(&"Sample_Name").unwrap().to_string());
        match record.get(&"Sample_Project") {
            Some(&project_name) if project_name.len() > 0 => {
                project_names.push(Some(project_name.to_string()))
            }
            Some(_) | None => project_names.push(None),
        }
        match record.get(&"Index") {
            Some(&idx) if idx.len() > 0 => sample_idx.push(idx.as_bytes().to_vec()),
            Some(_) | None => (),
        }
        match record.get(&"Index2") {
            Some(&idx2) if idx2.len() > 0 => sample_idx2.push(idx2.as_bytes().to_vec()),
            Some(_) | None => (),
        }
    }

    let sample_data: HashMap<_, _> = lanes
        .iter()
        .map(|(&i, (sample_names, project_names, idx_vec, idx2_vec))| {
            (
                i,
                make_sample_maps(sample_names, project_names, idx_vec, idx2_vec, max_distance),
            )
        })
        .collect();

    Ok(sample_data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use maplit::{hashmap, hashset};
    use std::fs;

    use ndarray::array;

    static ROOT: &str = "test_data/sample_data";

    #[test]
    fn with_conflict_no_index2_w_lanes() {
        let samplesheet = PathBuf::from(ROOT).join("w_conflict_no_index2_w_lanes.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let test_contents = fs::read_to_string("test_data/hamming_distance_1_test.txt").unwrap();

        let set0 = hashset! {0};

        let expected_lane1_index: HashMap<_, _> = test_contents
            .split_whitespace()
            .map(|ix| (ix.as_bytes().to_vec(), set0.clone()))
            .collect();

        let expected_lane1 = Samples {
            sample_names: vec!["sample_1".to_string()],
            project_names: vec![None],
            index_vec: vec![hashset! { vec![65, 67, 84, 71, 67, 71, 65, 65] }],
            index_map: expected_lane1_index,
            index2_vec: Vec::new(),
            index2_map: HashMap::new(),
        };

        let test_contents = fs::read_to_string("test_data/hamming_distance_1_test2.txt").unwrap();

        let expected_lane2_index: HashMap<_, _> = test_contents
            .split_whitespace()
            .map(|ix| (ix.as_bytes().to_vec(), set0.clone()))
            .collect();

        let expected_lane2 = Samples {
            sample_names: vec!["sample_2".to_string()],
            project_names: vec![None],
            index_vec: vec![hashset! { vec![65, 67, 84, 67, 65, 84, 67, 67] }],
            index_map: expected_lane2_index,
            index2_vec: Vec::new(),
            index2_map: HashMap::new(),
        };

        let mut expected_sampledata = HashMap::new();

        expected_sampledata.insert(1, expected_lane1);
        expected_sampledata.insert(2, expected_lane2);

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn with_conflict_w_index2_w_lanes() {
        let samplesheet = PathBuf::from(ROOT).join("w_conflict_w_index2_w_lanes.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let expected_index = hashmap! {
            b"CCCCT".to_vec() => hashset! {0},
            b"CCCCC".to_vec() => hashset! {1},
        };

        let expected_index2 = hashmap! {
            b"AAAAT".to_vec() => hashset! {0},
            b"AAAAA".to_vec() => hashset! {1},
        };

        let lane_map = sampledata.get(&1).unwrap();

        assert_eq!(lane_map.index_map, expected_index);
        assert_eq!(lane_map.index2_map, expected_index2);
    }

    #[test]
    fn with_conflict_in_separate_lanes_w_index2() {
        let samplesheet = PathBuf::from(ROOT).join("w_conflict_in_separate_lanes_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let set0 = hashset! {0};

        let expected_index: HashMap<_, _> = vec![
            "CCTCT", "CGCCT", "CCCCT", "CCCCA", "CCACT", "CNCCT", "GCCCT", "CCGCT", "ACCCT",
            "TCCCT", "CTCCT", "CCCGT", "NCCCT", "CCNCT", "CCCTT", "CCCCG", "CACCT", "CCCCN",
            "CCCNT", "CCCCC", "CCCAT",
        ]
        .iter()
        .map(|i| (i.as_bytes().to_vec(), set0.clone()))
        .collect();

        let expected_index2: HashMap<_, _> = vec![
            "AAAGA", "AGAAA", "AAATA", "AAAAG", "AACAA", "TAAAA", "AAANA", "CAAAA", "NAAAA",
            "AANAA", "ATAAA", "AAACA", "ACAAA", "AAAAN", "GAAAA", "AATAA", "AAGAA", "AAAAA",
            "AAAAT", "AAAAC", "ANAAA",
        ]
        .iter()
        .map(|i| (i.as_bytes().to_vec(), set0.clone()))
        .collect();

        assert_eq!(sampledata.len(), 2);
        for lane_map in sampledata.values() {
            assert_eq!(lane_map.index_map, expected_index);
            assert_eq!(lane_map.index2_map, expected_index2);
        }
    }

    #[test]
    fn no_conflict_w_index2() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let project_names = vec![Some("project_1".to_string()), Some("project_2".to_string())];

        let set0 = hashset! {0};
        let set1 = hashset! {1};

        let mut expected_index = HashMap::new();

        for ix in vec![
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC", "NGGGG", "GGGCG", "GGNGG",
            "GCGGG", "GNGGG", "GGGNG", "GGGGA", "GGGGN", "GAGGG", "GGGTG", "GGGGG", "GGCGG",
            "CGGGG", "GTGGG", "GGGGT",
        ] {
            expected_index.insert(ix.as_bytes().to_vec(), set0.clone());
        }

        for ix in vec![
            "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT", "TTTTT", "NTTTT", "ATTTT",
            "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT", "TTTTA", "TTTTC", "TTTTG", "CTTTT",
            "TTTGT", "TTTTN", "TTATT",
        ] {
            expected_index.insert(ix.as_bytes().to_vec(), set1.clone());
        }

        let mut expected_index2 = HashMap::new();

        for ix2 in vec![
            "ATAAA", "CAAAA", "AAGAA", "AACAA", "AAAAC", "AAAAG", "AAAAA", "AAATA", "ANAAA",
            "AAAGA", "AGAAA", "GAAAA", "NAAAA", "ACAAA", "AANAA", "AAAAN", "AAAAT", "AAACA",
            "AAANA", "TAAAA", "AATAA",
        ] {
            expected_index2.insert(ix2.as_bytes().to_vec(), set0.clone());
        }

        for ix2 in vec![
            "CCGCC", "CCNCC", "CCTCC", "TCCCC", "CGCCC", "CCCNC", "ACCCC", "CCCTC", "CCCGC",
            "CCCAC", "CCCCA", "CACCC", "GCCCC", "CCCCC", "CCCCT", "CTCCC", "CNCCC", "NCCCC",
            "CCCCG", "CCACC", "CCCCN",
        ] {
            expected_index2.insert(ix2.as_bytes().to_vec(), set1.clone());
        }

        let mut expected_sampledata = SampleData::new();

        expected_sampledata.insert(
            0,
            Samples {
                sample_names,
                project_names,
                index_map: expected_index,
                index_vec: vec![
                    hashset! { vec![71, 71, 71, 71, 71] },
                    hashset! { vec![84, 84, 84, 84, 84] },
                ],
                index2_map: expected_index2,
                index2_vec: vec![
                    hashset! { vec![65, 65, 65, 65, 65] },
                    hashset! { vec![67, 67, 67, 67, 67] },
                ],
            },
        );

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn double_index_samples() {
        let samplesheet = PathBuf::from(ROOT).join("double_index_samples.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();

        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let project_names = vec![Some("project_1".to_string()), Some("project_2".to_string())];

        let set0 = hashset! {0};
        let set1 = hashset! {1};

        let mut expected_index = HashMap::new();

        for ix in vec![
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC", "NGGGG", "GGGCG", "GGNGG",
            "GCGGG", "GNGGG", "GGGNG", "GGGGA", "GGGGN", "GAGGG", "GGGTG", "GGGGG", "GGCGG",
            "CGGGG", "GTGGG", "GGGGT", "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT",
            "TTTTT", "NTTTT", "ATTTT", "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT", "TTTTA",
            "TTTTC", "TTTTG", "CTTTT", "TTTGT", "TTTTN", "TTATT",
        ] {
            expected_index.insert(ix.as_bytes().to_vec(), set0.clone());
        }

        for ix2 in vec![
            "ATAAA", "CAAAA", "AAGAA", "AACAA", "AAAAC", "AAAAG", "AAAAA", "AAATA", "ANAAA",
            "AAAGA", "AGAAA", "GAAAA", "NAAAA", "ACAAA", "AANAA", "AAAAN", "AAAAT", "AAACA",
            "AAANA", "TAAAA", "AATAA", "CCGCC", "CCNCC", "CCTCC", "TCCCC", "CGCCC", "CCCNC",
            "ACCCC", "CCCTC", "CCCGC", "CCCAC", "CCCCA", "CACCC", "GCCCC", "CCCCC", "CCCCT",
            "CTCCC", "CNCCC", "NCCCC", "CCCCG", "CCACC", "CCCCN",
        ] {
            expected_index.insert(ix2.as_bytes().to_vec(), set1.clone());
        }

        let mut expected_sampledata = SampleData::new();

        expected_sampledata.insert(
            0,
            Samples {
                sample_names,
                project_names,
                index_map: expected_index,
                index_vec: vec![
                    hashset! { vec![71, 71, 71, 71, 71], vec![84, 84, 84, 84, 84] },
                    hashset! { vec![65, 65, 65, 65, 65], vec![67, 67, 67, 67, 67] },
                ],
                index2_map: HashMap::new(),
                index2_vec: Vec::new(),
            },
        );

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn sample_lookup() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();
        let lane = &sampledata.get(&0).unwrap();

        let idx1 = array![71, 71, 71, 71, 71];
        let idx1a = array![71, 84, 71, 71, 71];
        let idx2 = array![65, 65, 65, 65, 65];

        let idx3 = array![84, 84, 84, 84, 84];
        let idx3a = array![84, 67, 84, 84, 84];
        let idx4 = array![67, 67, 67, 67, 67];

        // good lookup, one index, exact match
        assert_eq!(Some((0, true)), lane.get_sample(&[idx1.view()]));

        // good lookup, one index, inexact match
        assert_eq!(Some((0, false)), lane.get_sample(&[idx1a.view()]));

        // good lookup, two indices, exact match
        assert_eq!(
            Some((1, true)),
            lane.get_sample(&[idx3.view(), idx4.view()])
        );

        // good lookup, two indices, exact match
        assert_eq!(
            Some((1, false)),
            lane.get_sample(&[idx3a.view(), idx4.view()])
        );

        // bad single index
        assert_eq!(None, lane.get_sample(&[idx2.view()]));

        // bad index, two indices
        assert_eq!(None, lane.get_sample(&[idx2.view(), idx1.view()]));

        // bad index2
        assert_eq!(None, lane.get_sample(&[idx1.view(), idx1.view()]));

        // indexes match different samples
        assert_eq!(None, lane.get_sample(&[idx1.view(), idx4.view()]));
    }

    #[test]
    fn any_sample_check() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();
        let lane = &sampledata.get(&0).unwrap();

        let idx1 = vec![71, 71, 71, 71, 71];
        let idx2 = vec![65, 65, 65, 65, 65];
        let idx3 = vec![84, 84, 84, 84, 84];
        let idx4 = vec![67, 67, 67, 67, 67];

        let idx1a = vec![71, 71, 71, 71, 65];
        let idx2g = vec![65, 65, 65, 65, 71];

        // any sample, one index
        assert!(lane.is_any_sample(&[idx1.clone()]));
        assert!(lane.is_any_sample(&[idx1a.clone()]));
        assert!(lane.is_any_sample(&[idx3.clone()]));

        // any sample, two indices
        assert!(lane.is_any_sample(&[idx1.clone(), idx2.clone()]));
        assert!(lane.is_any_sample(&[idx1a.clone(), idx2g.clone()]));
        assert!(lane.is_any_sample(&[idx3.clone(), idx4.clone()]));

        // no sample, one index
        assert!(!lane.is_any_sample(&[idx2.clone()]));

        // no sample, two indices
        assert!(!lane.is_any_sample(&[idx2.clone(), idx1.clone()]));
        assert!(!lane.is_any_sample(&[idx1a.clone(), idx3.clone()]));
        assert!(!lane.is_any_sample(&[idx4.clone(), idx4.clone()]));
    }

    #[test]
    fn make_sample_maps() {
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let project_names = vec![Some("project_1".to_string()), Some("project_2".to_string())];
        let mut expected_index = HashMap::new();

        let set0 = hashset! {0};
        let set1 = hashset! {1};

        for ix in vec![
            "TGGGG", "GGTGG", "AGGGG", "GGGAG", "GGAGG", "GGGGC", "NGGGG", "GGGCG", "GGNGG",
            "GCGGG", "GNGGG", "GGGNG", "GGGGA", "GGGGN", "GAGGG", "GGGTG", "GGGGG", "GGCGG",
            "CGGGG", "GTGGG", "GGGGT",
        ] {
            expected_index.insert(ix.as_bytes().to_vec(), set0.clone());
        }

        for ix in vec![
            "TCTTT", "TGTTT", "TTGTT", "TATTT", "TTCTT", "TTNTT", "TTTTT", "NTTTT", "ATTTT",
            "GTTTT", "TTTNT", "TNTTT", "TTTAT", "TTTCT", "TTTTA", "TTTTC", "TTTTG", "CTTTT",
            "TTTGT", "TTTTN", "TTATT",
        ] {
            expected_index.insert(ix.as_bytes().to_vec(), set1.clone());
        }

        let index_vec = vec![b"GGGGG".to_vec(), b"TTTTT".to_vec()];

        let actual_mapping =
            super::make_sample_maps(&sample_names, &project_names, &index_vec, &[], 1);

        assert_eq!(actual_mapping.index_map, expected_index);
    }

    #[test]
    fn make_sample_maps_conflict() {
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let project_names = vec![None, None];
        let index_vec = vec![b"ACTG".to_vec(), b"ACTC".to_vec()];

        let expected_index = hashmap! {
            b"ACTG".to_vec() => hashset! {0},
            b"ACTC".to_vec() => hashset! {1},
        };

        let actual_mapping =
            super::make_sample_maps(&sample_names, &project_names, &index_vec, &[], 1);

        assert_eq!(actual_mapping.index_map, expected_index);
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn no_file() {
        let samplesheet = PathBuf::from(ROOT).join("no_file.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"invalid UTF-8 in field 0 near byte index 0"#)]
    fn bad_file() {
        let samplesheet = PathBuf::from(ROOT).join("bad_file.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"No samples found in samplesheet"#)]
    fn empty_file() {
        let samplesheet = PathBuf::from("test_data/empty_file");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"Can't demux two different samples using the same indices"#)]
    fn sample_collision() {
        let samplesheet = PathBuf::from(ROOT).join("sample_collision.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"Samplesheet does not have a Sample_Name column"#)]
    fn no_sample() {
        let samplesheet = PathBuf::from(ROOT).join("no_sample_name.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"Samplesheet does not have an Index column"#)]
    fn no_index() {
        let samplesheet = PathBuf::from(ROOT).join("no_index.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"Samplesheet is missing index2 for some samples"#)]
    fn missing_index() {
        let samplesheet = PathBuf::from(ROOT).join("missing_index.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"Samplesheet is missing project names for some samples"#)]
    fn missing_project() {
        let samplesheet = PathBuf::from(ROOT).join("missing_project.csv");
        read_samplesheet(samplesheet, 1).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"Got 3 indices?!"#)]
    fn weird_sample_lookup() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();
        let lane = &sampledata.get(&0).unwrap();

        lane.get_sample(&[
            array![71, 84].view(),
            array![65, 65].view(),
            array![65, 65].view(),
        ]);
    }

    #[test]
    #[should_panic(expected = r#"Got 3 indices?!"#)]
    fn weird_any_sample_check() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1).unwrap();
        let lane = &sampledata.get(&0).unwrap();

        lane.is_any_sample(&[vec![71, 84], vec![65, 65], vec![65, 65]]);
    }
}
