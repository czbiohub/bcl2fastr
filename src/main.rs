#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_xml_rs;
extern crate glob;

use std::path::Path;

mod parser;
mod filter_decoder;
mod locs_decoder;
mod cbcl_header_decoder;
mod extract_reads;

fn main() {
    // let filename_info = Path::new("src/test_data/RunInfo.xml");
    // let filename_params = Path::new("src/test_data/test_runparams.xml")

    // let cbcl_file = Path::new("/usr/src/bcl2fastr_testlane/C117.1/L001_1.cbcl");
    // let locs_file = Path::new("src/test_data/test_locs.locs");
    // let filter_file = Path::new("src/test_data/test_filter.filter");
    // parser::parse_run_info(filename_info);
    // parser::parse_run_params(filename_params);
    // cbcl_header_decoder::cbcl_header_decoder(cbcl_file);
    // locs_decoder::locs_decoder(locs_file);
    // filter_decoder::filter_decoder(filter_file);

    let locs_path = Path::new("src/test_data/test_locs.locs");
    let run_info_path = Path::new("src/test_data/RunInfo.xml");
    let run_params_path = Path::new("src/test_data/test_runparams.xml");
    let lane_path = Path::new("/usr/src/bcl2fastr_testlane/");

    extract_reads(locs_path, run_info_path, run_params_path, lane_path);


}
