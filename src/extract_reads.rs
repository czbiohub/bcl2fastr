use std::{
    fs::read_dir,
    path::Path,
};
use glob::glob;

use parser::parse_run_info;
use parser::parse_run_params;
use filter_decoder::filter_decoder;
use locs_decoder::locs_decoder;
use cbcl_header_decoder::cbcl_header_decoder;


static LANE_PARTS : u8 = 2;  // supports 2 parts per lane

// extracts the reads for a particular lane and is then able to pass those into processes
pub fn extract_reads(locs_path: &Path, run_info_path: &Path, run_params_path: &Path, lane_path: &Path) {

    // read in metadata for the run: locs path, run info path, run params paths
    let locs = locs_decoder(locs_path);
    let run_info = parse_run_info(run_info_path);
    let run_params = parse_run_params(run_params_path);

    // expected number of cycles based on manual inputs into Run Parameters
    let num_cycles = run_params.read1_cycles + run_params.read2_cycles
        + run_params.index1_cycles + run_params.index2_cycles;
    println!("Expected number of cycles: {}", num_cycles);
    // read in metadata for the lane: cbcl headers, filters:
    // read in cbcl and filter paths using glob
    let cbcl_paths = glob(lane_path, "/C*").
            expect("Failed to read glob pattern for C* dirs");
    let filter_paths = glob(lane_path, "*.filter").
            expect("Failed to read glob pattern for *.filter files");

    println!("{:#}", cbcl_paths);
    println!("{:#}", filter_paths);

    // confirm that expected num_cycles matches the actual ones you're getting
    println!("Actual number of cycles: {}", cbcl_paths.len());
    if cbcl_paths.len() != num_cycles {
        panic!("Expected number of cycles, {0}, does not match the number of
            directories output by the sequencer, {1}.", num_cycles, cbcl_paths.len());
    }

    // read in all of the header information for the cbcl files in this lane
    let mut headers = [[CBCLHeader; LANE_PARTS]; num_cycles];
    for (c, cbcl_path) in cbcl_paths.enumerate() {
        let lane_part_paths = fs::read_dir(cbcl_path).unwrap();

        // check if number of lane part cbcls match the hard-coded lane parts
        if lane_part_paths.len() != LANE_PARTS {
            panic!("Number of lane parts for CBCL directory {0}, {1} does not match the number of
                expected lane parts, {2}", c, lane_part_paths.len(), LANE_PARTS);
        }

        for (p, part_path) in lane_part_paths.enumerate() {
            headers[c][p] = cbcl_header_decoder(part_path);
        }
    }
    println!("{:#?}", headers);

    // read in the filter paths as filter structs (there should be one for each tile)
    let mut filters = [Filter; filter_paths.len()];
    for (f, filter_path) in filter_paths.enumerate() {
        filters[f] = filter_decoder(filter_path);
    }
    println!("{:#?}", filters);



}
