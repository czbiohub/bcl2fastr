use std::{
    Collections::HashMap,
    fs::File,
    fs,
    path::Path,
};

mod parser;
mod cbcl_decoder;
mod filter_decoder;
mod locations;

pub struct Lane {
    pub cbcl_paths : Vec<String>, // paths to each cbcl file
    pub headers : Vec<CBCLHeader>, // use cbcl_fils to read in headers
    pub tiles : Vec<i32>, // tile numbers (all of the tiles in this lane, can be read from any cbcl file (should be same))
    pub filter_paths : Vec<String>, // paths to each filter_file
    pub filters : HashMap<i32, Filter>, // [tile_number, Filter object] read filters using filter_decoder
    pub base_matrix : Vec<Vec<u8>>, // (0 - 4)
    pub qscore_matrix : Vec<Vec<u8>>, // (0 - 4)
    pub reads : Vec<Read>, // vec of Read structs built using the base_matrix and qscore_matrix as inputs
    pub indices : Vec<Index>, // keep track of the indices when you first read it? could be helpful?

}

impl Lane {
    // lane_path directory ex: "<run_id>/Data/Intensities/BaseCalls/L001" for Lane 1
    // run_id ex: 190414_A00111_0296_AHJCWWDSXX

    // reads in paths from the lane_path dir and sorts into subdirectories (cbcl) or files (filters)
    fn read_subpaths(lane_path : String) -> (Vec<String>, Vec<String>) {
        let path = Path::new(lane_path);
        let cbcl_paths = Vec::new();
        let filter_paths = Vec::new();
        for entry in fs::read_dir(path) {
            let subpath = entry.path();
            let metadata = fs::metadata(&subpath)?;
            if metadata.is_file() {
                filter_paths.push(subpath.into_os_string().into_string().unwrap();)
            } else {
                cbcl_paths.push(subpath.into_os_string().into_string().unwrap();)
            }
        }
        return (cbcl_paths, filter_paths);
    }

    fn read_cbcl_headers(cbcl_files) -> Vec<CBCLHeader> {

    }

    // get from run_info
    fn read_tile_nums(lane_path) -> Vec<i32> {

    }

    fn read_filters(lane_path, tiles) -> HashMap<i32, Filter> {

    }

    fn decode_base() -> u8 {

    }

    fn decode_qscore() -> u8 {

    }

    // output base_matrix, qscore_matrix
    fn decode_cbcls(cbcl_files) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    // remember to check if you have an odd # of reads that pass filter (in that case the last half
    // byte of the file doesn't have any real base/qscore info)
    fn apply_filters(filters, base_matrix, q_score_matrix) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    fn collect_reads(filters, base_matrix, q_score_matrix) {

    }

}

pub fn lane_collect(lane_path: String) -> Lane {
    let f = File::open(lane_path).unwrap();
    let (cbcl_paths, filter_paths) = Lane::read_subpaths(lane_path);

    let lane = Lane {
        cbcl_paths,

    };

    println!("{:#?}", lane);
    return lane;
}
