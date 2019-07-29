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
    pub headers : Vec<CBCLHeader>, // use cbcl_files to read in headers
    pub tiles : Vec<Vec<Tile>>, // tiles (holds all of the bases from each of the tiles), read in by cbcl_file
    // Tile struct holds base_matrix and qscore_matrix from each tile (sorted by Lane and Lane-part)
    // make sure to also store the tile number inside this struct
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
    fn read_subpaths(lane_path : &Path) -> (Vec<Path>, Vec<Path>) {
        let cbcl_paths = Vec::new();
        let filter_paths = Vec::new();
        for entry in fs::read_dir(lane_path) {
            let subpath = entry.path();
            let metadata = fs::metadata(&subpath)?;
            if metadata.is_file() {
                filter_paths.push(subpath)
            } else {
                cbcl_paths.push(subpath)
            }
        }
        return (cbcl_paths, filter_paths);
    }

    fn read_cbcl_header(cbcl_path : &Path) -> CBCLHeader {
        return cbcl_decoder::cbcl_decoder(cbcl_path)
    }

    fn read_cbcl_headers(cbcl_paths : Vec<&Path>) -> Vec<CBCLHeader> {
        let mut headers = Vec::new();
        for cbcl_path in cbcl_paths {
            let header = read_cbcl_header(cbcl_path);
            headers.push(header);
        }
        return headers
    }

    // read in Tile struct from one cbcl_file (cbcl_path)
    fn read_tiles(cbcl_path : &Path) -> Tile {

    }

    fn read_filters(lane_path : &Path, tiles : Vec<Vec<Tile>>) -> HashMap<i32, Filter> {

    }

    fn decode_base() -> u8 {

    }

    fn decode_qscore() -> u8 {

    }

    // output base_matrix, qscore_matrix
    fn compile_tiles(tiles : Vec<Vec<Tile>>) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    // remember to check if you have an odd # of reads that pass filter (in that case the last half
    // byte of the file doesn't have any real base/qscore info)
    fn apply_filters(filters, base_matrix, q_score_matrix) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {

    }

    fn collect_reads(filters, base_matrix, q_score_matrix) {

    }

}

pub fn lane_collect(lane_path: &Path) -> Lane {
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
