use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
};

#[derive(Debug, PartialEq)]
pub struct Locs {
    pub locs : Vec<Vec<f32>>, //f32
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


pub fn locs_decoder(locs_path: String) -> Locs {
    let f = File::open(locs_path).unwrap();
    let locs = Locs::from_reader(f).unwrap();
    println!("{:#?}", locs);
    return locs;
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_locs() {
        let test_file = "src/test_data/test_locs.locs".to_string();
        let actual_locs : Locs = locs_decoder(test_file);
        let expected_locs =
            Locs {
                locs : vec![
                    vec![
                        0.0,
                        0.0,
                    ],
                    vec![
                        1.8080986,
                        0.0,
                    ],
                    vec![
                        3.616197,
                        0.0,
                    ],
                    vec![
                        5.4242954,
                        0.0,
                    ],
                    vec![
                        7.232394,
                        0.0,
                    ],
                ],
            };
        assert_eq!(actual_locs, expected_locs)
    }
}
