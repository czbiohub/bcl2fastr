//! Read `*.filter` files into vectors of boolean values.

use byteorder::{LittleEndian, ReadBytesExt};

use std::{fs::File, io::Read, path::Path};

/// A filter is a vector of bytes representing pairs of booleans,
/// e.g. (false, false) = 0, (true, false) = 2, etc
pub type Filter = Vec<u8>;

/// Decode a `.filter` file into a `Filter` struct or panic
///
/// Format of a `.filter` file:
///  1. Two `u32` containing header info (ignored)
///  2. `u32` representing the number of clusters
///  3. `[u8; num_clusters]` of true/false (1 or 0) values
pub fn filter_decoder(filter_path: &Path) -> std::io::Result<Filter> {
    let mut rdr = File::open(filter_path)?;

    let _ = rdr.read_u64::<LittleEndian>()?;

    let num_clusters = rdr.read_u32::<LittleEndian>()? as usize;

    let mut bin_mask = vec![0u8; num_clusters];
    rdr.read_exact(&mut bin_mask)?;
    // if length is off, add extra 0
    if num_clusters % 2 == 1 {
        bin_mask.push(0);
    }

    let filter = bin_mask.chunks_exact(2).map(|x| 2 * x[0] + x[1]).collect();

    Ok(filter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode() {
        let test_file = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX")
            .join("Data/Intensities/BaseCalls/L001/s_1_1101.filter");
        let actual_filter = filter_decoder(&test_file).unwrap();
        let expected_filter = vec![
            0, 1, 3, 3, 0, 2, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 1, 3, 3, 3,
            2, 3, 3, 1, 3, 3, 2, 3, 2, 1, 3, 2, 3, 3, 3, 2, 3, 0, 3, 2, 0,
        ];
        assert_eq!(actual_filter, expected_filter);
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn no_file() {
        let test_file = Path::new("test_data/no_file.filter");
        filter_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"failed to fill whole buffer"#)]
    fn empty_file() {
        let test_file = Path::new("test_data/empty_file");
        filter_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"failed to fill whole buffer"#)]
    fn bad_8_bytes() {
        let test_file = Path::new("test_data/bad_data_8.bin");
        filter_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"failed to fill whole buffer"#)]
    fn bad_12_bytes() {
        let test_file = Path::new("test_data/bad_data_12.bin");
        filter_decoder(test_file).unwrap();
    }
}
