//! Reads `*.locs` files into vectors of float arrays

use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};


#[derive(Debug, PartialEq)]
/// Each element is an array of [x, y] locations, one for each cluster in a tile
pub struct Locs {
    pub locs: Vec<[f32; 2]>,
}


impl Locs {
    /// Creates a `Locs` struct from a file reader
    /// 
    /// Format of a `.locs` file: 
    ///  1. Two `u32` containing header info (ignored)
    ///  2. `u32` representing number of clusters (locations)
    ///  3. `[f32; 2 * num_clusters]` of [x, y] pairs
    fn from_reader(mut rdr: impl Read) -> io::Result<Self> {
        let _ = rdr.read_u64::<LittleEndian>()?;
        let num_clusters = rdr.read_u32::<LittleEndian>()? as usize;

        let mut loc_buffer = vec![0f32; num_clusters * 2];
        rdr.read_f32_into::<LittleEndian>(&mut loc_buffer)?;

        let mut locs = Vec::new();
        for loc_chunk in loc_buffer.chunks_exact(2) {
            let mut loc = [0f32; 2];
            loc.clone_from_slice(loc_chunk);
            locs.push(loc);
        }

        Ok(Locs { locs })
    }
}


/// Decode a `.locs` file into a `Locs` struct or panic
pub fn locs_decoder(locs_path: &Path) -> Locs {
    let f = match File::open(locs_path) {
        Err(e) => panic!(
            "couldn't open {}: {}", locs_path.display(), e
        ),
        Ok(file) => file,
    };

    let locs = match Locs::from_reader(f) {
        Err(e) => panic!(
            "Error reading locs from {}: {}", locs_path.display(), e
        ),
        Ok(locs) => locs,
    };

    locs
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_locs() {
        let test_file = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs");
        let actual_locs = locs_decoder(test_file);
        let expected_locs =
            Locs {
                locs: vec![
                    [0.0, 0.0], [1.8080986, 0.0], [3.616197, 0.0], [5.4242954, 0.0],
                    [7.232394, 0.0], [9.040493, 0.0], [10.848592, 0.0], [12.656691, 0.0],
                    [14.464789, 0.0], [16.272888, 0.0], [18.080986, 0.0], [19.889084, 0.0],
                    [21.697182, 0.0], [23.50528, 0.0], [25.313377, 0.0], [27.121475, 0.0],
                    [28.929573, 0.0], [30.73767, 0.0], [32.54577, 0.0], [34.353867, 0.0],
                    [36.161964, 0.0], [37.970062, 0.0], [39.77816, 0.0], [41.586258, 0.0],
                    [43.394356, 0.0], [45.202454, 0.0], [47.01055, 0.0], [48.81865, 0.0],
                    [50.626747, 0.0], [52.434845, 0.0], [54.242943, 0.0], [56.05104, 0.0],
                    [57.85914, 0.0], [59.667236, 0.0], [61.475334, 0.0], [63.283432, 0.0],
                    [65.09153, 0.0], [66.89963, 0.0], [68.707726, 0.0], [70.51582, 0.0],
                    [72.32392, 0.0], [74.13202, 0.0], [75.94012, 0.0], [77.748215, 0.0],
                    [79.55631, 0.0], [81.36441, 0.0], [83.17251, 0.0], [84.980606, 0.0],
                    [86.788704, 0.0], [88.5968, 0.0], [90.4049, 0.0], [92.213, 0.0],
                    [94.021095, 0.0], [95.82919, 0.0], [97.63729, 0.0], [99.44539, 0.0],
                    [101.25349, 0.0], [103.061584, 0.0], [104.86968, 0.0], [106.67778, 0.0],
                    [108.48588, 0.0], [110.293976, 0.0], [112.10207, 0.0], [113.91017, 0.0],
                    [115.71827, 0.0], [117.52637, 0.0], [119.334465, 0.0], [121.14256, 0.0],
                    [122.95066, 0.0], [124.75876, 0.0], [126.56686, 0.0], [128.37495, 0.0],
                    [130.18306, 0.0], [131.99117, 0.0], [133.79927, 0.0], [135.60738, 0.0],
                    [137.41548, 0.0], [139.22359, 0.0], [141.0317, 0.0], [142.8398, 0.0],
                    [144.6479, 0.0], [146.45601, 0.0], [148.26411, 0.0], [150.07222, 0.0],
                    [151.88033, 0.0], [153.68843, 0.0], [155.49654, 0.0], [157.30464, 0.0],
                    [159.11275, 0.0], [160.92085, 0.0], [162.72896, 0.0], [164.53706, 0.0],
                    [166.34517, 0.0], [168.15327, 0.0], [169.96138, 0.0], [171.76949, 0.0],
                    [173.57759, 0.0], [175.3857, 0.0], [177.1938, 0.0], [179.0019, 0.0]
                ]
            };
        assert_eq!(actual_locs, expected_locs)
    }

    #[test]
    #[should_panic(
      expected = r#"couldn't open test_data/no_file.locs: No such file or directory (os error 2)"#
    )]
    fn test_locs_no_file() {
        let test_file = Path::new("test_data/no_file.locs");
        locs_decoder(test_file);
    }
}
