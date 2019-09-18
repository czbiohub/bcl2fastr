use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

#[derive(Debug, PartialEq)]
pub struct Locs {
    pub locs : Vec<Vec<f32>>, //f32, tuples of xy coordinates corresponding to each tile
}


impl Locs {
    fn from_reader(mut rdr: impl Read) -> io::Result<Self> {
        let _ = rdr.read_u64::<LittleEndian>()?;
        let num_clusters = rdr.read_u32::<LittleEndian>()?;
        let mut locs = Vec::new();
        for _c in 0..num_clusters {
            let x = rdr.read_f32::<LittleEndian>()?;
            let y = rdr.read_f32::<LittleEndian>()?;
            locs.push(vec![x, y]);
        }

        Ok(Locs {
            locs,
        })
    }
}


pub fn locs_decoder(locs_path: &Path) -> Locs {
    let f = File::open(locs_path).unwrap();
    let locs = Locs::from_reader(f).unwrap();
    return locs;
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_locs() {
        let test_file = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs");
        let actual_locs : Locs = locs_decoder(test_file);
        let expected_locs =
            Locs {
                locs: vec![
                    vec![0.0, 0.0],
                    vec![1.8080986, 0.0], vec![3.616197, 0.0], vec![5.4242954, 0.0],
                    vec![7.232394, 0.0], vec![9.040493, 0.0], vec![10.848592, 0.0],
                    vec![12.656691, 0.0], vec![14.464789, 0.0], vec![16.272888, 0.0],
                    vec![18.080986, 0.0], vec![19.889084, 0.0], vec![21.697182, 0.0],
                    vec![23.50528, 0.0], vec![25.313377, 0.0], vec![27.121475, 0.0],
                    vec![28.929573, 0.0], vec![30.73767, 0.0], vec![32.54577, 0.0],
                    vec![34.353867, 0.0], vec![36.161964, 0.0], vec![37.970062, 0.0],
                    vec![39.77816, 0.0], vec![41.586258, 0.0], vec![43.394356, 0.0],
                    vec![45.202454, 0.0], vec![47.01055, 0.0], vec![48.81865, 0.0],
                    vec![50.626747, 0.0], vec![52.434845, 0.0], vec![54.242943, 0.0],
                    vec![56.05104, 0.0], vec![57.85914, 0.0], vec![59.667236, 0.0],
                    vec![61.475334, 0.0], vec![63.283432, 0.0], vec![65.09153, 0.0],
                    vec![66.89963, 0.0], vec![68.707726, 0.0], vec![70.51582, 0.0],
                    vec![72.32392, 0.0], vec![74.13202, 0.0], vec![75.94012, 0.0],
                    vec![77.748215, 0.0], vec![79.55631, 0.0], vec![81.36441, 0.0],
                    vec![83.17251, 0.0], vec![84.980606, 0.0], vec![86.788704, 0.0],
                    vec![88.5968, 0.0], vec![90.4049, 0.0], vec![92.213, 0.0],
                    vec![94.021095, 0.0], vec![95.82919, 0.0], vec![97.63729, 0.0],
                    vec![99.44539, 0.0], vec![101.25349, 0.0], vec![103.061584, 0.0],
                    vec![104.86968, 0.0], vec![106.67778, 0.0], vec![108.48588, 0.0],
                    vec![110.293976, 0.0], vec![112.10207, 0.0], vec![113.91017, 0.0],
                    vec![115.71827, 0.0], vec![117.52637, 0.0], vec![119.334465, 0.0],
                    vec![121.14256, 0.0], vec![122.95066, 0.0], vec![124.75876, 0.0],
                    vec![126.56686, 0.0], vec![128.37495, 0.0], vec![130.18306, 0.0],
                    vec![131.99117, 0.0], vec![133.79927, 0.0], vec![135.60738, 0.0],
                    vec![137.41548, 0.0], vec![139.22359, 0.0], vec![141.0317, 0.0],
                    vec![142.8398, 0.0], vec![144.6479, 0.0], vec![146.45601, 0.0],
                    vec![148.26411, 0.0], vec![150.07222, 0.0], vec![151.88033, 0.0],
                    vec![153.68843, 0.0], vec![155.49654, 0.0], vec![157.30464, 0.0],
                    vec![159.11275, 0.0], vec![160.92085, 0.0], vec![162.72896, 0.0],
                    vec![164.53706, 0.0], vec![166.34517, 0.0], vec![168.15327, 0.0],
                    vec![169.96138, 0.0], vec![171.76949, 0.0], vec![173.57759, 0.0],
                    vec![175.3857, 0.0], vec![177.1938, 0.0], vec![179.0019, 0.0]
                ]
            };
        assert_eq!(actual_locs, expected_locs)
    }
}
