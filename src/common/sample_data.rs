//! Represents SampleSheet information as struct with lanes
//! that have index to sample mappings with all indices included within distance 1
//! of original index or distance 0 if overlapping indices are present within distance 1

use std::path::PathBuf;
use std::collections::HashMap;

use crate::hamming_set::{check_conflict, hamming_set, singleton_set};


/// SampleData maps from lane number to the index maps for the lane. The maps are
/// chunked into different pieces, each corresponding to a set of samples that will be
/// processed together
pub type SampleData = HashMap<usize, Vec<Samples>>;


/// Lane has one or two maps from index to a tuple of (sample name, corrected index)
/// If there is only one index, index2 will be empty and shouldn't be used
#[derive(Debug, PartialEq)]
pub struct Samples {
    pub sample_names: Vec<String>,
    index: HashMap<Vec<u8>, (String, String)>,
    index2: HashMap<Vec<u8>, (String, String)>
}


impl Samples {
    /// Look up a sample given a vector of indices (themselves vectors)
    pub fn get_sample(&self, indices: &[Vec<u8>]) -> Option<(String, String)> {
        match indices.len() {
            1 => self.get_1index_sample(&indices[0]),
            2 => self.get_2index_sample(&indices[0], &indices[1]),
            x => panic!("Got {} indices?!", x)
        }
    }

    /// helper function for when there is one index
    fn get_1index_sample(&self, i: &[u8]) -> Option<(String, String)> {
        match self.index.get(i) {
            Some(sample_t) => return Some(sample_t.clone()),
            None => return None
        };
    }

    /// helper function for when there are two indices
    fn get_2index_sample(&self, i: &[u8], i2: &[u8]) -> Option<(String, String)> {
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


/// Function to go from a lane worth of sample and index vectors to a Lane struct
/// which will include the necessary error-correction, up to some limit `max_distance`
fn make_sample_maps(
    sample_names: &[String],
    index_vec: &[Vec<u8>],
    index2_vec: &[Vec<u8>],
    max_distance: usize,
    n_threads: usize
) -> Vec<Samples> {
    // index_vec should be full
    assert_eq!(sample_names.len(), index_vec.len(), "Missing indexes for some samples");

    // index2_vec should either be full or empty, nothing in between
    assert!(
        index2_vec.len() == index_vec.len() || index2_vec.len() == 0,
        "Samplesheet is missing index2 for some samples"
    );

    // start at distance 0: just map samples to indices
    let mut index_hash_sets: Vec<_> = index_vec.iter().map(singleton_set).collect();
    let mut index2_hash_sets: Vec<_> = index2_vec.iter().map(singleton_set).collect();

    if check_conflict(&sample_names, &index_hash_sets, &index2_hash_sets) {
        panic!("Can't demux two different samples using the same indices");
    }

    for _ in 1..=max_distance {
        let new_index_hash_sets: Vec<_> = index_hash_sets.iter()
            .map(hamming_set)
            .collect();

        let new_index2_hash_sets: Vec<_> = index2_hash_sets.iter()
            .map(hamming_set)
            .collect();

        if check_conflict(&sample_names, &new_index_hash_sets, &new_index2_hash_sets) {
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

    // this is empty if there is no second index
    let index2_vec: Vec<_> = index2_vec.iter()
        .cloned()
        .map( |index| unsafe { String::from_utf8_unchecked(index) })
        .collect();

    // now that we have the final index set, we need to chunk the samples,
    // making sure to put all samples with the same name in the same chunk
    let mut samples = HashMap::new();
    let mut samples2 = HashMap::new();

    for ((sample_name, og_index), hash_set) in sample_names.iter()
        .zip(index_vec.iter())
        .zip(index_hash_sets) {

        let index_map = samples.entry(sample_name).or_insert(Vec::new());
        for index in hash_set.iter().cloned() {
            index_map.push((index, (sample_name.clone(), og_index.clone())));
        }
    }

    for ((sample_name, og_index2), hash_set) in sample_names.iter()
        .zip(index2_vec.iter())
        .zip(index2_hash_sets) {

        let index2_map = samples2.entry(sample_name).or_insert(Vec::new());
        for index2 in hash_set.iter().cloned() {
            index2_map.push((index2, (sample_name.clone(), og_index2.clone())));
        }
    }

    // distribute samples over chunks 
    let chunk_size = std::cmp::max(samples.len() / n_threads, 1);

    let lane_vec: Vec<_> = sample_names.chunks(chunk_size).map(
        |sample_chunk|
        {
            let index_map: HashMap<_, _> = sample_chunk.iter().filter_map(
                |sample_name| samples.get(sample_name)
            ).flatten().cloned().collect();

            let index2_map: HashMap<_, _> = sample_chunk.iter().filter_map(
                |sample_name| samples2.get(sample_name)
            ).flatten().cloned().collect();

            Samples { 
                sample_names: sample_chunk.to_vec(), 
                index: index_map, 
                index2: index2_map 
            }
        }
    ).collect();

    lane_vec
}


/// loads a sample sheet and converts it into a SampleData struct. Our version
/// automatically determines the mismatch rate that prevents conflicts, up to
/// a specified maximum
pub fn read_samplesheet(
    samplesheet: PathBuf, max_distance: usize, n_threads: usize
) -> std::io::Result<SampleData> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_path(samplesheet)?;

    // ignore any rows before the Data section
    let rows: Vec<_> = rdr.records()
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

    for record in rows[2..].iter()
        .map( |r| 
            rows[1].iter()
                .zip(r.iter())
                .collect::<HashMap<_, _>>()
        ) {
        let lane: usize = match record.get(&"Lane") {
            Some(lane) => lane.parse().unwrap(),
            None => 0,
        };

        let (sample_name, sample_idx, sample_idx2) = lanes.entry(lane).or_insert_with(
            || (Vec::new(), Vec::new(), Vec::new())
        );

        sample_name.push(record.get(&"Sample_Name").unwrap().to_string());
        match record.get(&"Index") {
            Some(&idx) if idx.len() > 0 => sample_idx.push(idx.as_bytes().to_vec()),
            Some(_) | None => ()
        }
        match record.get(&"Index2") {
            Some(&idx2) if idx2.len() > 0 => sample_idx2.push(idx2.as_bytes().to_vec()),
            Some(_) | None => ()
        };
    }

    let sample_data: HashMap<_, _> = lanes.iter()
        .map( |(&i, (smp_names, idx_vec, idx2_vec))|
            (i, make_sample_maps(smp_names, idx_vec, idx2_vec, max_distance, n_threads))
        )
        .collect();
 
    Ok(sample_data)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    static ROOT: &str = "test_data/sample_data";

    #[test]
    fn with_conflict_no_index2_w_lanes() {
        let samplesheet = PathBuf::from(ROOT).join("w_conflict_no_index2_w_lanes.csv");
        let sampledata = read_samplesheet(samplesheet, 1, 2).unwrap();

        let test_contents = fs::read_to_string(
            "test_data/hamming_distance_1_test.txt"
        ).unwrap();

        let expected_lane1_index: HashMap<_, _> = test_contents.split_whitespace()
            .map( |ix| 
                (
                    ix.as_bytes().to_vec(),
                    ("sample_1".to_string(), "ACTGCGAA".to_string())
                )
            )
            .collect();

        let expected_lane1 = Samples {
            sample_names: vec!["sample_1".to_string()],
            index: expected_lane1_index, 
            index2: HashMap::new()
        };

        let test_contents = fs::read_to_string(
            "test_data/hamming_distance_1_test2.txt"
        ).unwrap();

        let expected_lane2_index: HashMap<_, _> = test_contents.split_whitespace()
            .map( |ix|
                (
                    ix.as_bytes().to_vec(),
                    ("sample_2".to_string(), "ACTCATCC".to_string())
                )
            )
            .collect();

        let expected_lane2 = Samples {
            sample_names: vec!["sample_2".to_string()],
            index: expected_lane2_index,
            index2: HashMap::new()
        };

        let mut expected_sampledata = HashMap::new();

        expected_sampledata.insert(1, vec![expected_lane1]);
        expected_sampledata.insert(2, vec![expected_lane2]);

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn with_conflict_w_index2_w_lanes() {
        let samplesheet = PathBuf::from(ROOT).join("w_conflict_w_index2_w_lanes.csv");
        let sampledata = read_samplesheet(samplesheet, 1, 1).unwrap();

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

        assert_eq!(lane_map[0].index, expected_index);
        assert_eq!(lane_map[0].index2, expected_index2);
    }

    #[test]
    fn with_conflict_in_separate_lanes_w_index2() {
        let samplesheet = PathBuf::from(ROOT).join(
            "w_conflict_in_separate_lanes_w_index2.csv"
        );
        let sampledata = read_samplesheet(samplesheet, 1, 1).unwrap();

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
            assert_eq!(lane_map[0].index, expected_index[lane - 1]);
            assert_eq!(lane_map[0].index2, expected_index2[lane - 1]);
        }
    }

    #[test]
    fn no_conflict_w_index2() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1, 1).unwrap();

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
            vec![
                Samples {
                    sample_names,
                    index: expected_index,
                    index2: expected_index2
                }
            ]
        );

        assert_eq!(sampledata, expected_sampledata);
    }

    #[test]
    fn sample_lookup() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1, 1).unwrap();
        let lane = &sampledata.get(&0).unwrap()[0];

        // good lookup, one index
        assert_eq!(
            lane.get_sample(&[vec![71, 84, 71, 71, 71]]).unwrap(), 
            ("sample_1".to_string(), "GGGGG".to_string())
        );

        // good lookup, two indices
        assert_eq!(
            lane.get_sample(
                &[vec![71, 84, 71, 71, 71], vec![65, 65, 65, 65, 65]]
            ).unwrap(), 
            ("sample_1".to_string(), "GGGGG+AAAAA".to_string())
        );

        // bad single index
        assert_eq!(
            lane.get_sample(&[vec![65, 65, 65, 65, 65]]),
            None
        );

        // bad index, two indices
        assert_eq!(
            lane.get_sample(
                &[vec![65, 65, 65, 65, 65], vec![71, 84, 71, 71, 71]]
            ),
            None
        );

        // bad index2
        assert_eq!(
            lane.get_sample(
                &[vec![71, 84, 71, 71, 71], vec![71, 84, 71, 71, 71]]
            ),
            None
        );

        // indexes match different samples
        assert_eq!(
            lane.get_sample(
                &[vec![71, 84, 71, 71, 71], vec![67, 67, 67, 67, 71]]
            ),
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
            &sample_names, &index_vec, &[], 1, 1
        );

        assert_eq!(actual_mapping[0].index, expected_index);
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
            &sample_names, &index_vec, &[], 1, 1
        );

        assert_eq!(actual_mapping[0].index, expected_index);
    }

    #[test]
    #[should_panic(
        expected = r#"No such file or directory"#
    )]
    fn no_file() {
        let samplesheet = PathBuf::from(ROOT).join("no_file.csv");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"invalid UTF-8 in field 0 near byte index 0"#
    )]
    fn bad_file() {
        let samplesheet = PathBuf::from(ROOT).join("bad_file.csv");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"No samples found in samplesheet"#
    )]
    fn empty_file() {
        let samplesheet = PathBuf::from("test_data/empty_file");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Can't demux two different samples using the same indices"#
    )]
    fn sample_collision() {
        let samplesheet = PathBuf::from(ROOT).join("sample_collision.csv");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Samplesheet does not have a Sample_Name column"#
    )]
    fn no_sample() {
        let samplesheet = PathBuf::from(ROOT).join("no_sample_name.csv");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Samplesheet does not have an Index column"#
    )]
    fn no_index() {
        let samplesheet = PathBuf::from(ROOT).join("no_index.csv");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Samplesheet is missing index2 for some samples"#
    )]
    fn missing_index() {
        let samplesheet = PathBuf::from(ROOT).join("missing_index.csv");
        read_samplesheet(samplesheet, 1, 1).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Got 3 indices?!"#
    )]
    fn weird_sample_lookup() {
        let samplesheet = PathBuf::from(ROOT).join("no_conflict_w_index2.csv");
        let sampledata = read_samplesheet(samplesheet, 1, 1).unwrap();
        let lane = &sampledata.get(&0).unwrap()[0];

        lane.get_sample(&[vec![71, 84], vec![65, 65], vec![65, 65]]).unwrap();
    }
}
