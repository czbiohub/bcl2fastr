use byteorder::{LittleEndian, ReadBytesExt}
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

    // return base matrix and q_score matrix for each tile
    fn decode_tile() -> (Vec<Vec<u8>>, Vec<Vec<u8>>){

    }

}
