use std::{
    Collections::HashMap,
    fs::File,
    fs,
    path::Path,
};
use glob::glob

mod parser;
mod cbcl_decoder;
mod filter_decoder;
mod locations;

pub struct Lane {
    pub lane_path : &Path,  // path to the directory for the lane
    pub cbcl_paths : Result<Paths, PatternError>, // paths to each cbcl file
    pub cbcls : Vec<CBCL>,  // holds each decoded CBCL struct for this lane
    pub tiles : Vec<Vec<Tile>>, // tiles (holds all of the bases from each of the tiles), read in by cbcl_file
    // Tile struct holds base_matrix and qscore_matrix from each tile (sorted by Lane and Lane-part)
    // make sure to also store the tile number inside this struct
    pub filter_paths : Result<Paths, PatternError>, // paths to each filter_file
    pub filters : HashMap<i32, Filter>, // [tile_number, Filter object] read filters using filter_decoder
    pub base_matrix : Vec<Vec<u8>>, // (0 - 4)
    pub qscore_matrix : Vec<Vec<u8>>, // (0 - 4)
    pub reads : Vec<Read>, // vec of Read structs built using the base_matrix and qscore_matrix as inputs
    pub indices : Vec<Index>, // keep track of the indices when you first read it? could be helpful?
    pub locations: Vec<Vec<f32>>,  // tuples of xy coordinates corresponding to each tile
}

impl Lane {
    // lane_path directory ex: "<run_id>/Data/Intensities/BaseCalls/L001" for Lane 1
    // run_id ex: 190414_A00111_0296_AHJCWWDSXX

    // reads in paths from the lane_path dir and sorts into subdirectories (cbcl) or files (filters)
    fn read_subpaths(&mut self, lane_path : &Path) -> (Result<Paths, PatternError>, Result<Paths, PatternError>) {
        let cbcl_paths = glob(lane_path + "/C*");
        let filter_paths = glob(lane_path + "*.filter");
        self.cbcl_paths = cbcl_paths;
        self.filter_paths = filter_paths;
    }

    fn read_cbcls(&mut self, cbcl_paths : Vec<&Path>) -> Vec<CBCLHeader> {
        let mut cbcls = Vec::new();
        for cbcl_path in cbcl_paths {
            let cbcl = CBCL::decode_cbcl(cbcl_path).unwrap();
            cbcls.push(cbcl);
        }
        self.cbcls = cbcls;
    }

    fn read_filters(lane_path : &Path, tiles : Vec<Vec<Tile>>) -> HashMap<i32, Filter> {

    }

    // remember to check if you have an odd # of reads that pass filter (in that case the last half
    // byte of the file doesn't have any real base/qscore info)
    fn apply_filters(filters, base_matrix, q_score_matrix) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    // use the base_matrix and qscore_matrix from each tile to compile the fastq files
    fn collect_reads(filters, base_matrix, q_score_matrix) {

    }

    fn decode_lane(&mut self, lane_path : &Path) -> Lane {
        self.lane_path = lane_path;
        let f = File::open(lane_path).unwrap();
        Lane::read_subpaths(&mut self);

        let lane = Lane {
            self.lane_path,
            self.cbcl_paths,
            self.filter_paths,
        }

        println!("{:#?}", lane)
    }

}

pub fn lane_collect(lane_path: &Path) -> Lane {
    self.lane_path =
    let f = File::open(lane_path).unwrap();
    let (cbcl_paths, filter_paths) = Lane::read_subpaths(lane_path);

    let lane = Lane {
        cbcl_paths,
        filter_paths,
    };

    println!("{:#?}", lane);
    return lane;
}


#[cfg(test)]
mod tests {

    use super::*;

    // **need to fill in this test once a test file exists***
    #[test]
    fn test_subpaths() {
        let test_path = Path::new(); // lane path
        let (actual_cbcl_paths : Vec<Path>, actual_filter_paths : Vec<Path>) = Lane::read_subpaths(test_path);
        let expected_cbcl_paths = vec![Path::new(""), Path::new(""), Path::new()];
        let expected_filter_paths = vec![Path::new(""), Path::new(""), Path::new()];
    }
}
