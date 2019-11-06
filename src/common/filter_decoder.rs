//! Read `*.filter` files into vectors of boolean values.

use byteorder::{LittleEndian, ReadBytesExt};

use std::{
    fs::File,
    io::Read,
    path::Path,
};


/// A filter is a vector of booleans representing whether each cluster
/// in a tile has passed quality filtering.
pub type Filter = Vec<bool>;


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

    let filter = bin_mask.into_iter().map(|x| x != 0).collect();

    Ok(filter)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode() {
        let test_file = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001/s_1_1101.filter");
        let actual_filter = filter_decoder(test_file).unwrap();
        let expected_filter = vec![
            false, false, false, true, true, true, true, true, false, false,
            true, false, true, true, true, true, true, true, true, true,
            true, false, true, true, true, true, true, true, true, true,
            true, true, true, true, true, true, true, true, true, true,
            false, true, true, true, true, true, true, true, true, true,
            false, true, true, true, true, true, true, true, true, false,
            true, true, true, true, false, true, true, true, true, true,
            true, false, true, true, true, false, false, true, true, true,
            true, false, true, true, true, true, true, true, true, false,
            true, true, false, false, true, true, true, false, false, false,
        ];
        assert_eq!(actual_filter, expected_filter);
    }

    #[test]
    #[should_panic(
      expected = r#"No such file or directory"#
    )]
    fn no_file() {
        let test_file = Path::new("test_data/no_file.filter");
        filter_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"failed to fill whole buffer"#
    )]
    fn empty_file() {
        let test_file = Path::new("test_data/empty_file");
        filter_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"failed to fill whole buffer"#
    )]
    fn bad_8_bytes() {
        let test_file = Path::new("test_data/bad_data_8.bin");
        filter_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"failed to fill whole buffer"#
    )]
    fn bad_12_bytes() {
        let test_file = Path::new("test_data/bad_data_12.bin");
        filter_decoder(test_file).unwrap();
    }
}
