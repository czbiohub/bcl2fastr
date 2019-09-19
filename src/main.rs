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
    let locs_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs");
    let run_info_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/RunInfo.xml");
    let lane_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001");

    extract_reads::extract_reads(locs_path, run_info_path, lane_path, vec![0, 1]).unwrap();
}
