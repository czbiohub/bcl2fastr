use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};


#[derive(Debug, PartialEq)]
pub struct Filter {
    pub bin_mask : Vec<bool>, //binary mask of booleans (one value per cluster in a tile)
}

impl Filter {
    fn from_reader(mut rdr: impl Read) -> io::Result<Self> {
        let _ = rdr.read_u64::<LittleEndian>()?;
        let num_clusters = rdr.read_u32::<LittleEndian>()?;
        let mut bin_mask = Vec::new();
        for _c in 0..num_clusters {
            let val = rdr.read_u8()?;
            bin_mask.push(val != 0);
        }

        Ok(Filter {
            bin_mask,
        })
    }
}

pub fn filter_decoder(filter_path: &Path) -> Filter {
    let f = match File::open(filter_path) {
        Err(e) => panic!(
            "couldn't open {}: {}", filter_path.display(), e
        ),
        Ok(file) => file,
    };

    let filter = match Filter::from_reader(f) {
        Err(e) => panic!(
            "Error reading filter from {}: {}", filter_path.display(), e
        ),
        Ok(f) => f,
    };

    filter
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_filter() {
        let test_file = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001/s_1_1101.filter");
        let actual_filter: Filter = filter_decoder(test_file);
        let expected_filter =
            Filter {
                bin_mask: vec![
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
                ],
            };
        assert_eq!(actual_filter, expected_filter);
    }

    #[test]
    #[should_panic(
      expected = r#"couldn't open test_data/no_file.filter: No such file or directory (os error 2)"#
    )]
    fn test_filter_no_file() {
        let test_file = Path::new("test_data/no_file.filter");
        filter_decoder(test_file);
    }
}
