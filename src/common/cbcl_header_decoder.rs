//! Read the header from CBCL file and decode into a struct of useful information about
//! the file, to allow efficient tile extraction later.

use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
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
    pub bins: Vec<u8>,
    pub num_tile_records: u32,
    pub tiles: Vec<u32>,
    pub non_pf_clusters_excluded: bool,
    pub start_pos: Vec<u64>,
    pub uncompressed_size: Vec<usize>,
    pub compressed_size: Vec<usize>,
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
    ///     
    ///     Note: we only store the tile number, and compute chunked block sizes
    ///  9. `u8` flag for whether this file is only reads that pass quality filtering
    pub fn from_path(cbcl_path: &Path, tile_chunk: usize) -> std::io::Result<Self> {
        let mut rdr = File::open(cbcl_path)?;

        let version = rdr.read_u16::<LittleEndian>()?;
        let header_size = rdr.read_u32::<LittleEndian>()?;
        let bits_per_basecall = rdr.read_u8()?;
        let bits_per_qscore = rdr.read_u8()?;

        assert_eq!(bits_per_basecall, 2);
        assert_eq!(bits_per_qscore, 2);

        let number_of_bins = rdr.read_u32::<LittleEndian>()?;
        let mut bin_buffer = vec![0u32; (2 * number_of_bins) as usize];
        rdr.read_u32_into::<LittleEndian>(&mut bin_buffer)?;

        let bins = bin_buffer
            .chunks_exact(2)
            .map(|bc| bc[1].max(2) as u8 + 33)
            .collect();

        let num_tile_records = rdr.read_u32::<LittleEndian>()?;
        let mut tile_buffer = vec![0u32; (4 * num_tile_records) as usize];
        rdr.read_u32_into::<LittleEndian>(&mut tile_buffer)?;

        let tile_offsets: Vec<[u32; 4]> = tile_buffer
            .chunks_exact(4)
            .map(|tc| [tc[0], tc[1], tc[2], tc[3]])
            .collect();

        let non_pf_clusters_excluded = rdr.read_u8()? != 0;

        let tiles = tile_offsets.iter().map(|t| t[0]).collect();

        let start_pos = tile_offsets
            .iter()
            .scan(header_size, |pos, &t| {
                *pos += t[3];
                Some(*pos - t[3])
            })
            .step_by(tile_chunk)
            .map(|v| v as u64)
            .collect();

        let uncompressed_size: Vec<usize> = tile_offsets
            .chunks(tile_chunk)
            .map(|c| c.iter().map(|v| v[2]).sum::<u32>() as usize)
            .collect();

        let compressed_size: Vec<usize> = tile_offsets
            .chunks(tile_chunk)
            .map(|c| c.iter().map(|v| v[3]).sum::<u32>() as usize)
            .collect();

        Ok(CBCLHeader {
            cbcl_path: cbcl_path.to_path_buf(),
            version,
            header_size,
            bits_per_basecall,
            bits_per_qscore,
            number_of_bins,
            bins,
            num_tile_records,
            tiles,
            non_pf_clusters_excluded,
            start_pos,
            uncompressed_size,
            compressed_size,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode() {
        let cbcl_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX")
            .join("Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl");
        let actual_cbclheader = CBCLHeader::from_path(&cbcl_path, 2).unwrap();
        let expected_cbclheader = CBCLHeader {
            cbcl_path: cbcl_path.to_path_buf(),
            version: 1,
            header_size: 97,
            bits_per_basecall: 2,
            bits_per_qscore: 2,
            number_of_bins: 4,
            bins: vec![35, 44, 58, 70],
            num_tile_records: 3,
            tiles: vec![1101, 1102, 1103],
            non_pf_clusters_excluded: false,
            start_pos: vec![97, 243],
            uncompressed_size: vec![100, 50],
            compressed_size: vec![146, 73],
        };
        assert_eq!(actual_cbclheader, expected_cbclheader)
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn no_file() {
        let cbcl_path = Path::new("test_data/no_file.cbcl");
        CBCLHeader::from_path(cbcl_path, 2).unwrap();
    }

    #[test]
    #[should_panic(expected = r#"failed to fill whole buffer"#)]
    fn bad_file() {
        let cbcl_path = Path::new("test_data/bad_data_8.bin");
        CBCLHeader::from_path(cbcl_path, 2).unwrap();
    }
}
