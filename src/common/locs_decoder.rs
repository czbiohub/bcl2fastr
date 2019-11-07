//! Reads `*.locs` files into vectors of float arrays

use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    path::Path,
};


/// Each element is an array of [x, y] locations, one for each cluster in a tile
pub type Locs = Vec<[u32; 2]>;


/// Decode a `.locs` file into a `Locs` struct or panic
/// 
/// Format of a `.locs` file: 
///  1. Two `u32` containing header info (ignored)
///  2. `u32` representing number of clusters (locations)
///  3. `[f32; 2 * num_clusters]` of [x, y] pairs
/// 
/// To go from f32 to the integer coordinates bcl2fastq outputs, we use the conversion
/// round((v as f64) * 10. + 1000.) as u32
pub fn locs_decoder(locs_path: &Path) -> std::io::Result<Locs> {
    let mut rdr = File::open(locs_path)?;

    let _ = rdr.read_u64::<LittleEndian>()?;

    let num_clusters = rdr.read_u32::<LittleEndian>()? as usize;

    let mut loc_buffer = vec![0f32; num_clusters * 2];

    rdr.read_f32_into::<LittleEndian>(&mut loc_buffer)?;

    // do the bananas conversion to get the right coordinates
    let loc_buffer: Vec<_> = loc_buffer.iter()
        .cloned()
        .map(|v| ((v as f64) * 10. + 1000.).round())
        .map(|v| v as u32)
        .collect();

    let locs = loc_buffer.chunks_exact(2).map(|v| [v[0], v[1]]).collect();

    Ok(locs)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode() {
        let test_file = Path::new(
            "test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs"
        );
        let actual_locs = locs_decoder(test_file).unwrap();
        let expected_locs = vec![
            [1000, 1000], [1018, 1000], [1036, 1000], [1054, 1000], [1072, 1000],
            [1090, 1000], [1108, 1000], [1127, 1000], [1145, 1000], [1163, 1000],
            [1181, 1000], [1199, 1000], [1217, 1000], [1235, 1000], [1253, 1000],
            [1271, 1000], [1289, 1000], [1307, 1000], [1325, 1000], [1344, 1000],
            [1362, 1000], [1380, 1000], [1398, 1000], [1416, 1000], [1434, 1000],
            [1452, 1000], [1470, 1000], [1488, 1000], [1506, 1000], [1524, 1000],
            [1542, 1000], [1561, 1000], [1579, 1000], [1597, 1000], [1615, 1000],
            [1633, 1000], [1651, 1000], [1669, 1000], [1687, 1000], [1705, 1000],
            [1723, 1000], [1741, 1000], [1759, 1000], [1777, 1000], [1796, 1000],
            [1814, 1000], [1832, 1000], [1850, 1000], [1868, 1000], [1886, 1000],
            [1904, 1000], [1922, 1000], [1940, 1000], [1958, 1000], [1976, 1000],
            [1994, 1000], [2013, 1000], [2031, 1000], [2049, 1000], [2067, 1000],
            [2085, 1000], [2103, 1000], [2121, 1000], [2139, 1000], [2157, 1000],
            [2175, 1000], [2193, 1000], [2211, 1000], [2230, 1000], [2248, 1000],
            [2266, 1000], [2284, 1000], [2302, 1000], [2320, 1000], [2338, 1000],
            [2356, 1000], [2374, 1000], [2392, 1000], [2410, 1000], [2428, 1000],
            [2446, 1000], [2465, 1000], [2483, 1000], [2501, 1000], [2519, 1000],
            [2537, 1000], [2555, 1000], [2573, 1000], [2591, 1000], [2609, 1000],
            [2627, 1000], [2645, 1000], [2663, 1000], [2682, 1000], [2700, 1000],
            [2718, 1000], [2736, 1000], [2754, 1000], [2772, 1000], [2790, 1000]
        ];
        assert_eq!(actual_locs, expected_locs)
    }

    #[test]
    #[should_panic(
        expected = r#"No such file or directory"#
    )]
    fn no_file() {
        let test_file = Path::new("test_data/no_file.locs");
        locs_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"failed to fill whole buffer"#
    )]
    fn empty_file() {
        let test_file = Path::new("test_data/empty_file");
        locs_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"failed to fill whole buffer"#
    )]
    fn bad_8_bytes() {
        let test_file = Path::new("test_data/bad_data_8.bin");
        locs_decoder(test_file).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"failed to fill whole buffer"#
    )]
    fn bad_12_bytes() {
        let test_file = Path::new("test_data/bad_data_12.bin");
        locs_decoder(test_file).unwrap();
    }
}
