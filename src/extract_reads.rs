extern crate flate2;
extern crate glob;

use std::{
    fs::File,
    io::prelude::*,
    io::SeekFrom,
    path::Path,
};
use flate2::read::MultiGzDecoder;

use crate::run_info_parser;
use crate::filter_decoder;
use crate::locs_decoder;
use crate::cbcl_header_decoder::{
    cbcl_header_decoder,
    CBCLHeader,
};


fn extract_tiles(header: &CBCLHeader, tile_idces: (usize, usize)) -> std::io::Result<Vec<u8>> {
    // identify first and last tile index to later get start and end byte position
    let (first_idx, last_idx) = tile_idces;

    // calculate start byte position
    let start_pos = (header.header_size + header.tile_offsets[0..first_idx].iter().map(|v| v[3]).sum::<u32>()) as u64;

    // calculate end byte position and expected size of the decompressed tile(s)
    // index 2 is the uncompressed tile size
    let uncompressed_size = header.tile_offsets[first_idx..=last_idx].iter().map(|v| v[2]).sum::<u32>() as usize;
    // index 3 is the compressed tile size
    let compressed_size = header.tile_offsets[first_idx..=last_idx].iter().map(|v| v[3]).sum::<u32>() as usize;

    // open file and read whole file into a buffer
    let mut cbcl = File::open(&header.cbcl_path)?;
    cbcl.seek(SeekFrom::Start(start_pos))?;

    let mut read_buffer = vec![0u8; compressed_size];
    cbcl.read_exact(&mut read_buffer)?;

    // use MultiGzDecoder to uncompress the number of bytes summed over the offsets of all tile_idces
    let mut uncomp_bytes = vec![0u8; uncompressed_size];
    let mut gz = MultiGzDecoder::new(&read_buffer[..]);
    gz.read_exact(&mut uncomp_bytes)?;

    Ok(uncomp_bytes)
}


// cbcl_paths is of the shape (num_cycles, lane_parts)
fn extract_base_matrix(
    headers: &Vec<Vec<CBCLHeader>>, tile_idces: (usize, usize)
) -> std::io::Result<()> {
    // iterate through the cbcl_paths by cycle and then by lane part
    let num_cycles = headers.len();
    for c in 0..num_cycles {
        let num_parts = headers[c].len();
        for p in 0..num_parts {
            let _uncomp_bytes = extract_tiles(&headers[c][p], tile_idces);
        }
    }
    Ok(())
}


// extracts the reads for a particular lane and is then able to pass those into processes
pub fn extract_reads(
    locs_path: &Path, run_info_path: &Path, lane_path: &Path, tile_idces: (usize, usize)
) -> std::io::Result<()> {

    // read in metadata for the run: locs path, run info path, run params paths
    let _locs = locs_decoder::locs_decoder(locs_path);
    let _run_info = run_info_parser::parse_run_info(run_info_path);

    // read in metadata for the lane: cbcl headers, filters:
    // read in cbcl and filter paths using glob

    let mut headers = Vec::new();

    for c_path in glob::glob(lane_path.join("C*").to_str().unwrap()).expect(
        "Failed to read glob pattern for C* dirs"
    ).filter_map(Result::ok) {
        let mut lane_part_headers = Vec::new();

        for part_path in glob::glob(c_path.join("*cbcl").to_str().unwrap()).expect(
            "Failed to read glob pattern for CBCL files"
        ).filter_map(Result::ok) {
            lane_part_headers.push(cbcl_header_decoder(&part_path));
        }

        headers.push(lane_part_headers);
    }

    // read in the filter paths as filter structs (there should be one for each tile)
    let mut filters = Vec::new();

    for filter_path in glob::glob(lane_path.join("*.filter").to_str().unwrap()).expect(
        "Failed to read glob pattern for *.filter files"
    ).filter_map(Result::ok) {
        filters.push(filter_decoder::filter_decoder(&filter_path));
    }

    extract_base_matrix(&headers, tile_idces)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_tiles() {
        let cbcl_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl");
        let cbcl_header = cbcl_header_decoder(cbcl_path);
        let tile_idces = (0, 1);
        let expected_bytes = vec![
            212, 254, 220, 221, 166, 108, 217, 232, 236, 221,
            157, 216, 220, 220, 205, 222, 140, 212, 157, 254,
            199, 221, 237, 185, 252, 199, 237, 253, 253, 68,
            237, 205, 199, 199, 237, 109, 205, 79, 200, 220,
            76, 253, 204, 253, 95, 223, 238, 78, 79, 206,
            220, 152, 220, 157, 255, 196, 207, 207, 133, 78,
            236, 222, 205, 254, 237, 204, 198, 218, 236, 204,
            206, 204, 214, 207, 222, 204, 201, 221, 103, 207,
            204, 196, 204, 88, 216, 205, 222, 251, 253, 206,
            206, 237, 223, 220, 205, 76, 220, 205, 232, 220
        ];

        let uncomp_bytes = extract_tiles(&cbcl_header, tile_idces).unwrap();

        assert_eq!(uncomp_bytes, expected_bytes)
    }

    #[test]
    fn test_extract_reads() {
        let locs_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs");
        let run_info_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/RunInfo.xml");
        let lane_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001");
        let tile_idces = (0, 1);

        extract_reads(locs_path, run_info_path, lane_path, tile_idces).unwrap()
    }
}
