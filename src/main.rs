#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_xml_rs;
extern crate glob;

use std::path::Path;

mod parser;
mod cbcl_decoder;
mod filter_decoder;
mod locations;
mod tile_decoder;

fn main() {
    // let filename_info = Path::new("src/test_data/RunInfo.xml");
    // let filename_params = Path::new("src/test_data/test_runparams.xml")

    let cbcl_file = Path::new("/usr/src/bcl2fastr_testlane/C116.1/L001_1.cbcl");
    // let locs_file = Path::new("src/test_data/test_locs.locs");
    // let filter_file = Path::new("src/test_data/test_filter.filter");
    // parser::parse_run_info(filename_info);
    // parser::parse_run_params(filename_params);
    cbcl_decoder::CBCLHeader::decode_cbcl_header(cbcl_file);
    // locations::locs_decoder(locs_file);
    // filter_decoder::filter_decoder(filter_file);
}
