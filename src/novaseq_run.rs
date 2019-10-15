//! Represents a NovaSeq sequencing run as a struct 
//! that can be shared across threads

use std::path::PathBuf;
use std::collections::HashMap;

use crate::run_info_parser::{RunInfo, parse_run_info};
use crate::filter_decoder::{Filter, filter_decoder};
use crate::locs_decoder::{Locs, locs_decoder};
use crate::cbcl_header_decoder::{CBCLHeader, cbcl_header_decoder};


/// Represents a sequencing run, including a bunch of metadata
/// and the headers of all the CBCL files.
pub struct NovaSeqRun {
    /// the root path of the sequencing run
    pub run_path: PathBuf,
    /// RunInfo object, stores the contents of RunInfo.xml
    pub run_info: RunInfo,
    /// a single universal locs array, same for every tile
    pub locs: Locs,
    /// a map from (lane, surface) tuples to vectors of filters
    pub filters: HashMap<(u32, u32), Vec<Filter>>,
    /// a map from (lane, surface) tuples to vectors of CBCL headers
    pub headers: HashMap<(u32, u32), Vec<CBCLHeader>>,
}


impl NovaSeqRun {
    /// Loads a NovaSeqRun from the given run path
    /// 
    /// TODO: Multithread this, lots of parallel tasks here!
    pub fn read_path(run_path: PathBuf) -> std::io::Result<NovaSeqRun> {
        let run_info = parse_run_info(&run_path.join("RunInfo.xml"))?;
        let locs = locs_decoder(&run_path.join("Data/Intensities/s.locs"))?;

        let mut headers = HashMap::new();
        let mut filters = HashMap::new();

        let n_lanes = run_info.runs.flow_cell_layout.lane_count;
        let n_surfaces = run_info.runs.flow_cell_layout.surface_count;
        let n_cycles = run_info.runs.reads.read.iter().map(|v| v.num_cycles).sum::<u64>();

        for lane in 1..=n_lanes {
            for surface in 1..=n_surfaces {
                let mut lane_surface_headers = Vec::new();
                let mut lane_surface_filters = Vec::new();

                for cycle in 1..=n_cycles {
                    let cbcl_path = run_path.join(
                        format!(
                            "Data/Intensities/BaseCalls/L{:03}/C{}.1/L{:03}_{}.cbcl",
                            lane,
                            cycle,
                            lane,
                            surface,
                        )
                    );
                    lane_surface_headers.push(cbcl_header_decoder(&cbcl_path)?);
                }

                // tile numbers are not stored by surface in RunInfo, so we are
                // taking advantage of the headers having the right names
                for tile_no in lane_surface_headers[0].tile_offsets.iter().map(|v| v[0]) {
                    let filter_path = run_path.join(
                        format!(
                            "Data/Intensities/BaseCalls/L{:03}/s_{}_{}.filter",
                            lane,
                            lane,
                            tile_no,
                        )
                    );
                    lane_surface_filters.push(filter_decoder(&filter_path)?);
                }

                headers.insert((lane, surface), lane_surface_headers);
                filters.insert((lane, surface), lane_surface_filters);
            }
        }

        let novaseq_run = NovaSeqRun {
            run_path: run_path,
            run_info: run_info,
            locs: locs,
            filters: filters,
            headers: headers,
        };

        Ok(novaseq_run)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let _novaseq_run = NovaSeqRun::read_path(run_path).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"No such file or directory"#
    )]
    fn no_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXXX");
        let _novaseq_run = NovaSeqRun::read_path(run_path).unwrap();
    }
}
