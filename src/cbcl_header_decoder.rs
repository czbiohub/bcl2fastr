//! Read the header from CBCL file and decode into a struct of useful information about
//! the file, to allow efficient tile extraction later.

use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};


#[derive(Debug, PartialEq)]
/// Represents the header information from a CBCL file
pub struct CBCLHeader {
    pub cbcl_path: PathBuf,
    pub version: u16,
    pub header_size: u32,
    pub bits_per_basecall: u8,
    pub bits_per_qscore: u8,
    pub number_of_bins: u32,
    pub bins: Vec<[u32; 2]>,
    pub num_tile_records: u32,
    // [tile number, num clusters, uncompressed block size, compressed block size]
    pub tile_offsets: Vec<[u32; 4]>,
    pub non_pf_clusters_excluded: bool,
}


impl CBCLHeader {
    /// Reads the beginning of a CBCL file and stores the header information
    /// 
    /// Structure of a CBCL header:
    ///  1. `u16` for file version
    ///  2. `u32` for header size, in bytes
    ///  3. `u8` for number of bits per basecall (usually 2)
    ///  4. `u8` for number of bits per quality score (usually 2)
    ///  5. `u32` for the number of quality score bins
    ///  6. `number_of_bins` pairs of `u32` key-value pairs for quality score bins.
    ///     First value is key, second is lower bound on the PHRED score.
    ///  7. `u32` for number of tiles in the file
    ///  8. `num_tile_records` arrays of 4 `u32` values:
    ///     1. Tile number
    ///     2. Number of clusters in the tile
    ///     3. Uncompressed block size
    ///     4. Compressed block size
    ///  9. `u8` flag for whether this file is only reads that pass quality filtering
    pub fn from_reader(cbcl_path: &Path, mut rdr: impl Read) -> io::Result<Self> {
        let version = rdr.read_u16::<LittleEndian>()?;
        let header_size = rdr.read_u32::<LittleEndian>()?;
        let bits_per_basecall = rdr.read_u8()?;
        let bits_per_qscore = rdr.read_u8()?;

        let number_of_bins = rdr.read_u32::<LittleEndian>()?;
        let mut bin_buffer = vec![0u32; (2 * number_of_bins) as usize];
        rdr.read_u32_into::<LittleEndian>(&mut bin_buffer)?;

        let mut bins = Vec::new();
        for bin_chunk in bin_buffer.chunks_exact(2) {
            let mut bin = [0u32; 2];
            bin.clone_from_slice(bin_chunk);
            bins.push(bin);
        }

        let num_tile_records = rdr.read_u32::<LittleEndian>()?;
        let mut tile_buffer = vec![0u32; (4 * num_tile_records) as usize];
        rdr.read_u32_into::<LittleEndian>(&mut tile_buffer)?;

        let mut tile_offsets = Vec::new();
        for tile_chunk in tile_buffer.chunks_exact(4) {
            let mut tile = [0u32; 4];
            tile.clone_from_slice(tile_chunk);
            tile_offsets.push(tile);
        }
        let non_pf_clusters_excluded = rdr.read_u8()? != 0;


        Ok(CBCLHeader {
            cbcl_path: cbcl_path.to_path_buf(),
            version,
            header_size,
            bits_per_basecall,
            bits_per_qscore,
            number_of_bins,
            bins,
            num_tile_records,
            tile_offsets,
            non_pf_clusters_excluded,
        })
    }
}


/// Decode a `.cbcl` header into a `CBCLHeader` struct or panic
pub fn cbcl_header_decoder(cbcl_path: &Path) -> CBCLHeader {
    let f = match File::open(cbcl_path) {
        Err(e) => panic!(
            "couldn't open {}: {}", cbcl_path.display(), e
        ),
        Ok(file) => file,
    };

    let cbcl_header = match CBCLHeader::from_reader(cbcl_path, f) {
        Err(e) => panic!(
            "Error reading header from {}: {}", cbcl_path.display(), e
        ),
        Ok(ch) => ch,
    };

    cbcl_header
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cbclheader() {
        let cbcl_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl");
        let actual_cbclheader = cbcl_header_decoder(cbcl_path);
        let expected_cbclheader =
            CBCLHeader {
                cbcl_path: cbcl_path.to_path_buf(),
                version: 1,
                header_size: 97,
                bits_per_basecall: 2,
                bits_per_qscore: 2,
                number_of_bins: 4,
                bins: vec![[0, 0], [1, 11], [2, 25], [3, 37]],
                num_tile_records: 3,
                tile_offsets: vec![[1101, 100, 50, 73], [1102, 100, 50, 73], [1103, 100, 50, 73]],
                non_pf_clusters_excluded: false,
            };
        assert_eq!(actual_cbclheader, expected_cbclheader)
    }

    #[test]
    #[should_panic(
      expected = r#"couldn't open test_data/no_file.cbcl: No such file or directory (os error 2)"#
    )]
    fn test_cbclheader_no_file() {
        let cbcl_path = Path::new("test_data/no_file.cbcl");
        cbcl_header_decoder(cbcl_path);
    }
}
