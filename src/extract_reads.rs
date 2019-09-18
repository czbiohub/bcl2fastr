extern crate flate2;

use std::{
    fs::File,
    io::prelude::*,
    io::SeekFrom,
    path::Path,
    path::PathBuf,
};
use crate::glob;
use flate2::read::MultiGzDecoder;

use crate::parser;
use crate::filter_decoder;
use crate::locs_decoder;
use crate::cbcl_header_decoder::{
    cbcl_header_decoder,
    CBCLHeader,
};


// cbcl_paths is of the shape (num_cycles, lane_parts)
fn extract_base_matrix(
    headers : &Vec<Vec<CBCLHeader>>, cbcl_paths : &Vec<Vec<PathBuf>>, tile_idces : Vec<u32>
) -> std::io::Result<()> {
    // initialize base matrix and qscore matrix (Vec::new())
    // let mut base_matrix : Vec<Vec<u32>> = Vec::new();
    // let mut qscore_matrix : Vec<Vec<u32>> = Vec::new();

    // iterate through the cbcl_paths by cycle and then by lane part
    let num_cycles = cbcl_paths.len();
    for c in 0..num_cycles {
        let num_parts = cbcl_paths[c].len();
        for p in 0..num_parts {
            // assign header and path
            let header: &CBCLHeader = &headers[c][p];
            let cbcl_path: &PathBuf = &cbcl_paths[c][p];
            
            // identify first and last tile index to later get start and end byte position
            let first_idx = tile_idces[0] as usize;
            let last_idx = (tile_idces[tile_idces.len() - 1] + 1) as usize;

            // calculate start byte position
            // if the first index is 0, edge case where you can't slice from 0..0,
            // so need to have a separate condition
            let start_pos;
            if first_idx == 0 {
                start_pos = header.header_size as u64;
            } else {
                start_pos = (header.header_size + header.tile_offsets[..first_idx].iter().map(|v| v[3]).sum::<u32>()) as u64;
            }

            // calculate end byte position and expected size of the decompressed tile(s)
            // if number of tiles is 1, then slicing returns a u32 instead of a vector slice,
            // so .iter() needs to be avoided
            let compressed_size;
            let uncompressed_size;

            if tile_idces.len() == 1 {
                // index 2 is the uncompressed tile size
                uncompressed_size = header.tile_offsets[first_idx][2] as usize;
                // index 3 is the compressed tile size
                compressed_size = header.tile_offsets[first_idx][3] as usize;
            } else {
                // index 2 is the uncompressed tile size
                uncompressed_size = header.tile_offsets[first_idx..last_idx].iter().map(|v| v[2]).sum::<u32>() as usize;
                // index 3 is the compressed tile size
                compressed_size = header.tile_offsets[first_idx..last_idx].iter().map(|v| v[3]).sum::<u32>() as usize;
            }

            // open file and read whole file into a buffer
            let mut cbcl = File::open(cbcl_path)?;
            cbcl.seek(SeekFrom::Start(start_pos))?;

            let mut read_buffer = vec![0u8; compressed_size];
            cbcl.read_exact(&mut read_buffer)?;

            // use MultiGzDecoder to uncompress the number of bytes summed over the offsets of all tile_idces
            let mut uncomp_bytes = vec![0u8; uncompressed_size];
            let mut gz = MultiGzDecoder::new(&read_buffer[..]);
            gz.read_exact(&mut uncomp_bytes)?;

            // check that size of decompressed tiles matches the size expected
            let actual_size = uncomp_bytes.len();
            if  actual_size != uncompressed_size {
                panic!("Decompressed tile(s) were expected to be {0} bytes long but were {1} bytes long", uncompressed_size, actual_size);
            }
        }
    }
    Ok(())
}

// dcbcl is decompressed matrix of bases that is passed in to a process
// this is the function implemented by a specific process (assigned a specific set of indices)
//fn get_read(dcbcl : Vec<Vec<u32>>, indices : (Vec<u32>, Vec<u32>), tile_idces : Vec<u32>, loc : u32) {
//}


// extracts the reads for a particular lane and is then able to pass those into processes
pub fn extract_reads(locs_path: &Path, run_info_path: &Path, lane_path: &Path, tile_idces : Vec<u32>) -> std::io::Result<()> {

    // read in metadata for the run: locs path, run info path, run params paths
    let _locs = locs_decoder::locs_decoder(locs_path);
    let _run_info = parser::parse_run_info(run_info_path);

    // read in metadata for the lane: cbcl headers, filters:
    // read in cbcl and filter paths using glob

    let mut cbcl_paths = Vec::new(); 
    let mut headers = Vec::new();
    let c_paths : Vec<_> = glob::glob(lane_path.join("C*").to_str().unwrap()).
            expect("Failed to read glob pattern for C* dirs").collect();

    for c_path in &c_paths {
        let lane_part_paths : Vec<_> = glob::glob(c_path.as_ref().unwrap().join("*").to_str().unwrap()).unwrap().collect();

        let mut lane_part_headers = Vec::new();
        let mut formatted_part_paths = Vec::new();
        for part_path in lane_part_paths {
            // part_path.as_ref() changes the value to be unwrapped from Option<T> to Option<&T> which doesn't assume ownership
            lane_part_headers.push(cbcl_header_decoder(part_path.as_ref().unwrap()));
            formatted_part_paths.push(part_path.unwrap());
        }

        headers.push(lane_part_headers);
        cbcl_paths.push(formatted_part_paths);
    }

    let filter_paths : Vec<_> = glob::glob(lane_path.join("*.filter").to_str().unwrap()).
            expect("Failed to read glob pattern for *.filter files").collect();

    // read in the filter paths as filter structs (there should be one for each tile)
    let mut filters = Vec::new();
    for filter_path in &filter_paths {
        filters.push(filter_decoder::filter_decoder(filter_path.as_ref().unwrap()));
    }

    extract_base_matrix(&headers, &cbcl_paths, tile_idces)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_reads() {
        let locs_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs");
        let run_info_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/RunInfo.xml");
        let lane_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001");
        let tile_idces = vec![0, 1];

        extract_reads(locs_path, run_info_path, lane_path, tile_idces).unwrap()
    }
}
