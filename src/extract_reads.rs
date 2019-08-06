extern crate flate2;

use std::{
    fs::File,
    io::{
        prelude::*,
        SeekFrom,
    },
    path::Path,
};
use crate::glob;
use flate2::read::GzDecoder;

use crate::parser;
use crate::filter_decoder;
use crate::locs_decoder;
use crate::cbcl_header_decoder::{
    cbcl_header_decoder,
    CBCLHeader,
};


static LANE_PARTS : u32 = 2;  // supports 2 parts per lane


// cbcl_paths is of the shape (num_cycles, lane_parts)
fn extract_base_matrix(headers : Vec<Vec<CBCLHeader>>, cbcl_paths : Vec<Vec<&Path>>, tile_idces : Vec<u32>) {
    // initialize base matrix and qscore matrix (Vec::new())
    let base_matrix = Vec::new();
    let qscore_matrix = Vec::new();

    // iterate through the cbcl_paths by cycle and then by lane part
    let num_cycles = cbcl_paths.len();
    for c in 0..num_cycles {
        let num_parts = cbcl_paths[c].len();
        for p in 0..num_parts {
            // assign header and path
            let header : CBCLHeader = headers[c][p];
            let cbcl_path : &Path = cbcl_paths[c][p];
            
            // open file and seek to (header_size + start of tile_offsets of all the previous tiles before the first tile_idx) to skip those bytes of the file
            let mut cbcl = File::open(cbcl_path).unwrap();
            let first_idx = tile_idces[0] as usize;
            // Note: must be Iterator to implement .sum()
            let start_pos = header.header_size + header.tile_offsets[0..first_idx][3].iter().sum();
            // Note: .into() converts header.header_size from u32 into u64 which is required by SeekFrom::Start()
            cbcl.seek(SeekFrom::Start(header.header_size.into()));

            // use GzDecoder to decompress the number of bytes summed over the offsets of all tile_idces
            let last_idx = tile_idces[tile_idces.len() - 1] as usize;
            let end_pos = header.tile_offsets[first_idx..last_idx][3].iter().sum();
            cbcl.seek(SeekFrom::End(end_pos));
            let mut tile_bytes = GzDecoder::new(cbcl).unwrap();

            // check that size of decompressed tiles matches the size expected
            let expected_size = header.tile_offsets[first_idx..last_idx][2].iter().sum() as usize;
            let actual_size = tile_bytes.stream_len();
            if  actual_size != expected_size {
                panic!("Decompressed tile(s) were expected to be {0} bytes long but were {1} bytes long", expected_size, actual_size);
            }


        }
    }
}

// dcbcl is decompressed matrix of bases that is passed in to a process
// this is the function implemented by a specific process (assigned a specific set of indices)
fn get_read(dcbcl : Vec<Vec<u32>>, indices : (Vec<u32>, Vec<u32>), tile_idces : Vec<u32>, loc : u32) {
    

}


// extracts the reads for a particular lane and is then able to pass those into processes
pub fn extract_reads(locs_path: &Path, run_info_path: &Path, run_params_path: &Path, lane_path: &Path, tile_idces : Vec<u32>) {

    // read in metadata for the run: locs path, run info path, run params paths
    let locs = locs_decoder::locs_decoder(locs_path);
    let run_info = parser::parse_run_info(run_info_path);
    let run_params = parser::parse_run_params(run_params_path);

    // expected number of cycles based on manual inputs into Run Parameters
    let num_cycles = run_params.read1_cycles + run_params.read2_cycles
        + run_params.index1_cycles + run_params.index2_cycles;

    // read in metadata for the lane: cbcl headers, filters:
    // read in cbcl and filter paths using glob

    let mut cbcl_paths = Vec::new(); 
    let mut headers = Vec::new();
    let c_paths : Vec<_> = glob::glob(lane_path.join("C*").to_str().unwrap()).
            expect("Failed to read glob pattern for C* dirs").collect();

    for c_path in &c_paths {
        let lane_part_paths : Vec<_> = glob::glob(c_path.as_ref().unwrap().join("*").to_str().unwrap()).unwrap().collect();
        // println!("Lane part paths: {:#?}", lane_part_paths);

        // // check if number of lane part cbcls match the hard-coded lane parts
        // // `as usize` converts a u32 into a usize type which allows for inequality comparison
        // if lane_part_paths.len() != LANE_PARTS as usize {
        //     panic!("Number of lane parts, {0}, for CBCL directory {1:#?}, does not match the number of
        //         expected lane parts, {2}", lane_part_paths.len(), c_path.as_ref().unwrap(), LANE_PARTS);
        // }

        let mut lane_part_headers = Vec::new();
        // need to borrow lane_part_paths as a reference in order to avoid using the value after move later in the code
        for part_path in &lane_part_paths {
            // &part_path.as_ref() changes the value to be unwrapped from &Option<T> to &Option<&T> which doesn't assume ownership
            // println!("{:#?}", &part_path.as_ref().unwrap());
            lane_part_headers.push(cbcl_header_decoder(&part_path.as_ref().unwrap()));
        }

        cbcl_paths.push(vec!(lane_part_paths));
        headers.push(lane_part_headers);
    }

    let filter_paths : Vec<_> = glob::glob(lane_path.join("*.filter").to_str().unwrap()).
            expect("Failed to read glob pattern for *.filter files").collect();
    
    // // confirm that expected num_cycles matches the actual ones you're getting
    // if headers.len() != num_cycles as usize {
    //     panic!("Expected number of cycles, {0}, does not match the number of
    //         directories present in the input files, {1}", num_cycles, headers.len());
    // }

    // read in the filter paths as filter structs (there should be one for each tile)
    let mut filters = Vec::new();
    for filter_path in &filter_paths {
        filters.push(filter_decoder::filter_decoder(&filter_path.as_ref().unwrap()));
    }
    println!("{:#?}", filters);

    // 
    // Note: could just get rid of headers matrix and get each header as you go through that cbcl file
    let flat_headers = headers.iter().flatten();
    let flat_paths = cbcl_paths.iter().flatten();
    for cbcl_data in &flat_headers.zip(&flat_paths) {
        let header = cbcl_data.0.as_ref().unwrap();
        let path = cbcl_data.1.as_ref().unwrap();
    
        // for a particular process:
        extract_base_matrix(header, path, tile_idces)
        get_read()
    } 
    


}
