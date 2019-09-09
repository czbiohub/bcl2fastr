use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};


#[derive(Debug, PartialEq)]
pub struct CBCLHeader {
    pub version : u16, // H
    pub header_size : u32, // I
    pub bits_per_basecall : u8, // B
    pub bits_per_qscore : u8, // B
    pub number_of_bins : u32, //I
    pub bins : Vec<Vec<u32>>, //I
    pub num_tile_records : u32, //I
    pub tile_offsets : Vec<Vec<u32>>, //I [tile number, num clusters, uncompressed
                                                // block size, compressed block size]
    pub non_PF_clusters_excluded : u8, //B, converted from u8 to bool

}


impl CBCLHeader {
    pub fn from_reader(mut rdr: impl Read) -> io::Result<Self> {
        let version = rdr.read_u16::<LittleEndian>()?;
        let header_size = rdr.read_u32::<LittleEndian>()?;
        let bits_per_basecall = rdr.read_u8()?;
        let bits_per_qscore = rdr.read_u8()?;
        let number_of_bins = rdr.read_u32::<LittleEndian>()?;
        let mut bins = Vec::new();
        for _b in 0..number_of_bins {
            let from = rdr.read_u32::<LittleEndian>()?;
            let to = rdr.read_u32::<LittleEndian>()?;
            bins.push(vec![from, to]);
        }
        let num_tile_records = rdr.read_u32::<LittleEndian>()?;
        let mut tile_offsets = Vec::new();
        for _t in 0..num_tile_records {
            let tile_number = rdr.read_u32::<LittleEndian>()?;
            let num_clusters = rdr.read_u32::<LittleEndian>()?;
            let uncomp_block_size = rdr.read_u32::<LittleEndian>()?;
            let comp_block_size = rdr.read_u32::<LittleEndian>()?;
            tile_offsets.push(vec![tile_number, num_clusters, uncomp_block_size, comp_block_size]);
        }
        let non_PF_clusters_excluded = rdr.read_u8()?;


        Ok(CBCLHeader {
            version,
            header_size,
            bits_per_basecall,
            bits_per_qscore,
            number_of_bins,
            bins,
            num_tile_records,
            tile_offsets,
            non_PF_clusters_excluded,
        })
    }

}


pub fn cbcl_header_decoder(cbcl_path : &Path) -> CBCLHeader {
    let f = File::open(cbcl_path).unwrap();
    let cbcl_header = CBCLHeader::from_reader(f).unwrap();
    println!("{:#?}", cbcl_header);
    return cbcl_header;
}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_cbclheader() {
        let test_file = Path::new("test_data/test_lane/C116.1/L001_1.cbcl");
        let actual_cbclheader : CBCLHeader = cbcl_header_decoder(test_file);
        let expected_cbclheader =
            CBCLHeader {
                version: 1,
                header_size: 65,
                bits_per_basecall: 2,
                bits_per_qscore: 2,
                number_of_bins: 4,
                bins: vec![
                    vec![
                        0,
                        0,
                    ],
                    vec![
                        1,
                        11,
                    ],
                    vec![
                        2,
                        25,
                    ],
                    vec![
                        3,
                        37,
                    ],
                ],
                num_tile_records: 1,
                tile_offsets: vec![
                    vec![
                        1101,
                        3366129,
                        1683065,
                        1088959,
                    ]
                ],
                non_PF_clusters_excluded: 1,
            };
        assert_eq!(actual_cbclheader, expected_cbclheader)
    }
}
