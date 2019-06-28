use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
};


#[derive(Debug, PartialEq)]
pub struct CBCLHeader {
    pub version : u16, // H
    pub header_size : u16, // I  
    pub bits_per_basecall : u16, // B   
    pub bits_per_qscore : u16, // B 
    pub number_of_bins : u16, //I
}


impl CBCLHeader {
    fn from_reader(mut rdr: impl Read) -> io::Result<Self> {
        let version = rdr.read_u16::<LittleEndian>()?;
        let header_size = rdr.read_u16::<LittleEndian>()?;
        let bits_per_basecall = rdr.read_u16::<LittleEndian>()?;
        let bits_per_qscore = rdr.read_u16::<LittleEndian>()?;
        let number_of_bins = rdr.read_u16::<LittleEndian>()?;        

        Ok(CBCLHeader {
            version,
            header_size,
            bits_per_basecall,
            bits_per_qscore,
            number_of_bins,
        }) 
    }

}


pub fn cbcl_decoder(cbcl_path: String) {

    let f = File::open(cbcl_path).unwrap();
    let cbcl = CBCLHeader::from_reader(f);
    println!("{:#?}", cbcl);

}
