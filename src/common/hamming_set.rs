//! Functions for building hamming sets (all strings within a given distance of a seed)
//! and checking for conflicts between them

use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use rayon::prelude::*;

/// given a set of indices, return a new set of indices which include everything
/// within hamming distance 1
pub fn hamming_set(index_set: &HashSet<Vec<u8>>) -> HashSet<Vec<u8>> {
    const NUCLEOTIDES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

    let new_set: HashSet<_> = index_set
        .iter()
        .cloned()
        .map(|index| {
            let mut this_set = HashSet::new();
            for i in 0..index.len() {
                for c in NUCLEOTIDES.iter().cloned() {
                    let mut new_index = index.clone();
                    new_index[i] = c;
                    this_set.insert(new_index);
                }
            }
            this_set
        })
        .flatten()
        .collect();

    new_set
}

/// Function to check for overlaps between the sets of sample indices. If there are
/// two indices, then an overlap between one is allowed as long as the second index
/// is sufficient to distinguish them.
pub fn check_conflict(
    index_sets: &HashMap<String, HashSet<Vec<u8>>>,
    index2_sets: &HashMap<String, HashSet<Vec<u8>>>,
) -> bool {
    let sample_clash: HashSet<_> = index_sets
        .iter()
        .tuple_combinations()
        .par_bridge()
        .filter_map(|((s1, hset1), (s2, hset2))| {
            if s1 != s2 && hset1.intersection(hset2).count() > 0 {
                Some((std::cmp::min(s1, s2), std::cmp::max(s1, s2)))
            } else {
                None
            }
        })
        .collect();

    if index2_sets.len() == 0 {
        return sample_clash.len() > 0;
    }

    let sample_clash2: HashSet<_> = index2_sets
        .iter()
        .tuple_combinations()
        .par_bridge()
        .filter_map(|((s1, hset1), (s2, hset2))| {
            if s1 != s2 && hset1.intersection(hset2).count() > 0 {
                Some((std::cmp::min(s1, s2), std::cmp::max(s1, s2)))
            } else {
                None
            }
        })
        .collect();

    return sample_clash.intersection(&sample_clash2).count() > 0;
}

#[cfg(test)]
mod test {
    use super::*;
    use maplit::{hashmap, hashset};
    use std::fs;

    #[test]
    fn singleton_set() {
        let index = b"ACTGCGAA".to_vec();
        let index_set = hashset! { index.clone() };

        assert_eq!(index_set.len(), 1);
        assert_eq!(index_set.get(&index), Some(&index));
    }

    #[test]
    fn hamming_set_distance_1() {
        let initial_set = hashset! { b"ACTGCGAA".to_vec() };
        let actual_hammingset = hamming_set(&initial_set);

        let test_contents = fs::read_to_string("test_data/hamming_distance_1_test.txt").unwrap();

        let expected_hammingset: HashSet<_> = test_contents
            .split_whitespace()
            .map(|i| i.as_bytes().to_vec())
            .collect();

        assert_eq!(actual_hammingset, expected_hammingset);
    }

    #[test]
    fn hamming_set_distance_2() {
        let initial_set = hashset! { b"ACTGCGAA".to_vec() };
        let actual_hammingset = hamming_set(&hamming_set(&initial_set));

        let test_contents = fs::read_to_string("test_data/hamming_distance_2_test.txt").unwrap();

        let expected_hammingset: HashSet<_> = test_contents
            .split_whitespace()
            .map(|i| i.as_bytes().to_vec())
            .collect();

        assert_eq!(actual_hammingset, expected_hammingset);
    }

    #[test]
    fn hamming_conflicts() {
        let index1 = hashset! { b"ACTGCGAA".to_vec() };
        let index2 = hashset! { b"ACTGCGAT".to_vec() };
        let index3 = hashset! { b"ACTGCCTT".to_vec() };

        let hammingset1 = hashmap! {
            "sample_1".to_string() => hamming_set(&index1),
            "sample_2".to_string() => hamming_set(&index2),
        };
        let hammingset2 = hashmap! {
            "sample_1".to_string() => hamming_set(&index2),
            "sample_2".to_string() => hamming_set(&index3),
        };
        let hammingset3 = hashmap! {
            "sample_1".to_string() => hamming_set(&index1),
            "sample_2".to_string() => hamming_set(&index3),
        };
        let empty_hammingset = HashMap::new();

        assert!(check_conflict(&hammingset1, &empty_hammingset));
        assert!(check_conflict(&hammingset2, &empty_hammingset));
        assert!(!check_conflict(&hammingset3, &empty_hammingset));
    }
}
