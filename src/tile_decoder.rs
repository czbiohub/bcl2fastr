use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};


#[derive(Debug, PartialEq)]
pub struct Tile {
    pub base_matrix : Vec<Vec<u8>>,
    pub qscore_matrix : Vec<Vec<u8>>,
}

impl Tile {
    // 
    // // byte to base (0-4)
    // fn decode_base() -> u8 {
    //
    // }
    //
    // // byte to qscore (0-4)
    // fn decode_qscore() -> u8 {
    //
    // }
    //
    // // return base matrix and q_score matrix for each tile
    // fn decode_tile() -> (Vec<Vec<u8>>, Vec<Vec<u8>>){
    //
    // }

}
