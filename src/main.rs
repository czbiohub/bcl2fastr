use std::path::Path;

mod run_info_parser;
mod filter_decoder;
mod locs_decoder;
mod cbcl_header_decoder;
mod extract_reads;

fn main() {
    let locs_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/s.locs");
    let run_info_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/RunInfo.xml");
    let lane_path = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/Data/Intensities/BaseCalls/L001");
    let tile_idces = (0, 1);

    extract_reads::extract_reads(locs_path, run_info_path, lane_path, tile_idces).unwrap();
}
