#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_xml_rs;

mod parser;
mod cbcl_decoder;
mod locations;

fn main() {
    let filename_info = "src/test_data/RunInfo.xml".to_string();
    let filename_params = "src/test_data/test_runparams.xml".to_string();
    let cbcl_file = "src/test_data/L001_1_cbcl_header.cbcl".to_string();
    let locs_file = "src/test_data/test_locs.locs".to_string();

    parser::parse_run_info(filename_info);
    parser::parse_run_params(filename_params);
    cbcl_decoder::cbcl_decoder(cbcl_file);
    locations::locs_decoder(locs_file);
}
