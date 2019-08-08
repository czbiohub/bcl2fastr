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
use flate2::read::{
    GzDecoder,
    MultiGzDecoder,
};

use crate::parser;
use crate::filter_decoder;
use crate::locs_decoder;
use crate::cbcl_header_decoder::{
    cbcl_header_decoder,
    CBCLHeader,
};


static LANE_PARTS : u32 = 2;  // supports 2 parts per lane


// cbcl_paths is of the shape (num_cycles, lane_parts)
fn extract_base_matrix(headers : &Vec<Vec<CBCLHeader>>, cbcl_paths : Vec<Vec<std::path::PathBuf>>, tile_idces : Vec<u32>) -> std::io::Result<()> {
    // initialize base matrix and qscore matrix (Vec::new())
    let mut base_matrix : Vec<Vec<u32>> = Vec::new();
    let mut qscore_matrix : Vec<Vec<u32>> = Vec::new();

    // iterate through the cbcl_paths by cycle and then by lane part
    let num_cycles = cbcl_paths.len();
    for c in 0..num_cycles {
        let num_parts = cbcl_paths[c].len();
        for p in 0..num_parts {
            // assign header and path
            let header = &headers[c][p];
            let cbcl_path = &cbcl_paths[c][p];
            
            // open file and read whole file into a buffer 
            let mut cbcl = File::open(cbcl_path)?;
            let mut whole_buffer = Vec::new();
            cbcl.read_to_end(&mut whole_buffer)?;

            // identify first and last tile index to later get start and end byte position
            let first_idx = tile_idces[0] as usize;
            let last_idx = tile_idces[tile_idces.len() - 1] as usize;

            // calculate start byte position
            // if the first index is 0, edge case where you can't slice from 0..0, so need to have a separate condition
            let start_pos;
            if first_idx == 0 {
                start_pos = 0;
            } else {
                start_pos = (header.header_size + &header.tile_offsets[..first_idx][3].iter().sum::<u32>()) as usize;
            }
            
            // calculate end byte position and expected size of the decompressed tile(s)
            // if number of tiles is 1, then slicing returns a u32 instead of a vector slice, so .iter() needs to be avoided
            let num_tiles = tile_idces.len();
            let end_pos;
            let expected_size;
            if num_tiles == 1 {
                end_pos = header.tile_offsets[last_idx][3] as usize; // index 3 is the compressed tile size 
                expected_size = header.tile_offsets[last_idx][2] as usize; // index 2 is the uncompressed tile size
            } else {
                end_pos = header.tile_offsets[first_idx..last_idx][3].iter().sum::<u32>() as usize; // index 3 is the compressed tile size
                expected_size = header.tile_offsets[first_idx..last_idx][2].iter().sum::<u32>() as usize; // index 2 is the compressed tile size
            }
            
            // slice buffer by appropriate start_pos and end_pos
            let tile_bytes = &whole_buffer[start_pos..end_pos];
            // println!("{}", end_pos - start_pos);
            // println!("{}", tile_bytes.len());

            // use GzDecoder to decompress the number of bytes summed over the offsets of all tile_idces
            let mut uncomp_bytes = Vec::new();
            let mut gz = MultiGzDecoder::new(tile_bytes);
            println!("built new decoder");
            gz.read_to_end(&mut uncomp_bytes)?;
            
            println!("finished decoding");

            // check that size of decompressed tiles matches the size expected
            let actual_size = uncomp_bytes.len();
            println!("{0} {1}", actual_size, expected_size);
            if  actual_size != expected_size {
                panic!("Decompressed tile(s) were expected to be {0} bytes long but were {1} bytes long", expected_size, actual_size);
            }
        }
    }
    Ok(())
}

// dcbcl is decompressed matrix of bases that is passed in to a process
// this is the function implemented by a specific process (assigned a specific set of indices)
fn get_read(dcbcl : Vec<Vec<u32>>, indices : (Vec<u32>, Vec<u32>), tile_idces : Vec<u32>, loc : u32) {
    

}


// extracts the reads for a particular lane and is then able to pass those into processes
pub fn extract_reads(locs_path: &Path, run_info_path: &Path, run_params_path: &Path, lane_path: &Path, tile_idces : Vec<u32>) -> std::io::Result<()> {

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
        let mut formatted_part_paths = Vec::new();
        // need to borrow lane_part_paths as a reference in order to avoid using the value after move later in the code
        for part_path in lane_part_paths {
            // &part_path.as_ref() changes the value to be unwrapped from &Option<T> to &Option<&T> which doesn't assume ownership
            // println!("{:#?}", &part_path.as_ref().unwrap());
            // let test = part_path.as_ref().unwrap();
            lane_part_headers.push(cbcl_header_decoder(part_path.as_ref().unwrap()));
            formatted_part_paths.push(part_path.unwrap());
        }

        headers.push(lane_part_headers);
        cbcl_paths.push(formatted_part_paths);
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
        filters.push(filter_decoder::filter_decoder(filter_path.as_ref().unwrap()));
    }
    println!("{:#?}", filters);

    // 
    // let flat_headers = headers.iter().flatten();
    // let flat_paths = cbcl_paths.iter().flatten();
    // for cbcl_data in &flat_headers.zip(&flat_paths) {
    //     let header = cbcl_data.0.as_ref().unwrap();
    //     let path = cbcl_data.1.as_ref().unwrap();
    
    //     // for a particular process:
    //     extract_base_matrix(header, path, tile_idces)
    //     get_read()
    // } 
    
    extract_base_matrix(&headers, cbcl_paths, tile_idces)

}
