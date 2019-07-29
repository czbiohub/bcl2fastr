#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_xml_rs;

use std::path::Path;

mod parser;
mod cbcl_decoder;
mod filter_decoder;
mod locations;

fn main() {
    let filename_info = Path::new("src/test_data/RunInfo.xml");
    let filename_params = Path::new("src/test_data/test_runparams.xml");
    let cbcl_file = Path::new("src/test_data/L001_1_cbcl_header.cbcl");
    let locs_file = Path::new("src/test_data/test_locs.locs");
    let filter_file = Path::new("src/test_data/test_filter.filter"sh);
    parser::parse_run_info(filename_info);
    parser::parse_run_params(filename_params);
    cbcl_decoder::cbcl_decoder(cbcl_file);
    locations::locs_decoder(locs_file);
    filter_decoder::filter_decoder(filter_file);
}
