use std::Collections::HashMap;

mod parser;
mod cbcl_decoder;
mod filter_decoder;
mod locations;

pub struct Lane {
    pub cbcl_files : Vec<&str>, // paths to each cbcl file
    pub headers : Vec<CBCLHeader>, // use cbcl_fils to read in headers
    pub tiles : Vec<i32>, // tile numbers (all of the tiles in this lane, can be read from any cbcl file (should be same))
    pub filters : HashMap<i32, Filter>, // [tile_number, Filter object] read filters using filter_decoder
    pub base_matrix : Vec<Vec<u8>>, // (0 - 4)
    pub qscore_matrix : Vec<Vec<u8>>, // (0 - 4)
    pub reads : Vec<Read>, // vec of Read structs built using the base_matrix and qscore_matrix as inputs
    pub indices : Vec<Index>, // keep track of the indices when you first read it? could be helpful?

}

impl Lane {
    fn read_cbcl_paths(lane_path) -> Vec<&str> {

    }

    fn read_cbcl_headers(cbcl_files) -> Vec<CBCLHeader> {

    }

    fn read_tile_nums(lane_path) -> Vec<i32> {

    }

    fn read_filters(lane_path, tiles) -> HashMap<i32, Filter> {

    }

    fn decode_base() -> u8 {

    }

    fn decode_qscore() -> u8 {

    }

    fn decode_cbcls(cbcl_files) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    // remember to check if you have an odd # of reads that pass filter (in that case the last half
    // byte of the file doesn't have any real base/qscore info)
    fn apply_filters(filters, base_matrix, q_score_matrix) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    fn collect_reads(filters, base_matrix, q_score_marix) {

    }


}
