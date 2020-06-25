//! Represents a NovaSeq sequencing run as a struct
//! that can be shared across threads

use std::collections::HashMap;
use std::path::PathBuf;

use log::{debug, info};
use rayon::prelude::*;

use crate::cbcl_header_decoder::CBCLHeader;
use crate::filter_decoder::{filter_decoder, Filter};
use crate::locs_decoder::{locs_decoder, Locs};
use crate::run_info_parser::{parse_run_info, RunInfo};

/// Represents a sequencing run, including a bunch of metadata
/// and the headers of all the CBCL files.
pub struct NovaSeqRun {
    /// the root path of the sequencing run
    pub run_path: PathBuf,
    /// RunInfo object, stores the contents of RunInfo.xml
    pub run_info: RunInfo,
    /// a string with the run info formatted for read headers
    pub run_id: String,
    /// a single universal locs array, same for every tile
    pub locs: Locs,
    /// a map from [lane, surface] tuples to vectors of filters
    pub filters: HashMap<[usize; 2], Vec<Filter>>,
    /// a map from [lane, surface] tuples to vectors of post-filter "filters"
    /// which are needed to remove empty half-bytes at the end of tiles
    pub pf_filters: HashMap<[usize; 2], Vec<Filter>>,
    /// a map from [lane, surface] to vectors of tiles, so we can produce
    /// the read header
    pub tile_ids: HashMap<[usize; 2], Vec<u32>>,
    /// a map from [lane, surface] to vectors of number of reads that pass filter,
    /// because we need this value a lot
    pub n_pfs: HashMap<[usize; 2], Vec<usize>>,
    /// a map from [lane, surface] to vectors of CBCL headers for the reads
    pub read_headers: HashMap<[usize; 2], Vec<Vec<CBCLHeader>>>,
    /// a map from [lane, surface] to vectors of CBCL headers for the indices
    pub index_headers: HashMap<[usize; 2], Vec<Vec<CBCLHeader>>>,
}

impl NovaSeqRun {
    /// Loads a NovaSeqRun located in `run_path`, aggregating filters and headers
    /// in chunks of `tile_chunk` tiles each. If `index_only` is true, will only
    /// load in data for cycles that are in indexes, and adjusts `index_ix` attribute
    /// accordingly. Uses threads to load the data in parallel.
    pub fn read_path(run_path: PathBuf, index_only: bool) -> std::io::Result<NovaSeqRun> {
        let run_info = parse_run_info(&run_path.join("RunInfo.xml"))?;
        let run_id = format!(
            "@{}:{}:{}",
            run_info.instrument, run_info.number, run_info.flowcell,
        );

        // need to repeat locs for each tile in tile_chunk
        let locs = locs_decoder(&run_path.join("Data/Intensities/s.locs"))?;

        let mut read_headers = HashMap::new();
        let mut index_headers = HashMap::new();

        let mut filters = HashMap::new();
        let mut pf_filters = HashMap::new();
        let mut tile_ids = HashMap::new();
        let mut n_pfs = HashMap::new();

        for lane in 1..=run_info.flowcell_layout.lane_count {
            for surface in run_info.flowcell_layout.surface_range.clone() {
                info!("lane {} - surface {}", lane, surface);

                let mut lane_surface_read_headers = Vec::new();
                let mut lane_surface_index_headers = Vec::new();

                for read in run_info.reads.iter() {
                    if index_only && !read.is_indexed_read {
                        continue;
                    }

                    let these_headers: Vec<CBCLHeader> = (read.start..read.end)
                        .into_par_iter()
                        .map(|cycle| {
                            let cbcl_path = run_path.join("Data/Intensities/BaseCalls").join(
                                format!("L{:03}/C{}.1/L{:03}_{}.cbcl", lane, cycle, lane, surface),
                            );

                            match CBCLHeader::from_path(&cbcl_path) {
                                Ok(header) => header,
                                Err(e) => {
                                    panic!("Error reading header {} {}", cbcl_path.display(), e)
                                }
                            }
                        })
                        .collect();

                    if read.is_indexed_read {
                        lane_surface_index_headers.push(these_headers);
                    } else {
                        lane_surface_read_headers.push(these_headers);
                    }
                }

                let mut lane_surface_filters = Vec::new();
                let mut lane_surface_tile_ids = Vec::new();

                // tile numbers are not stored by surface in RunInfo, so we are
                // taking advantage of the headers having the right names.
                // We will always have index_headers, but not always read_headers
                lane_surface_index_headers[0][0]
                    .tiles
                    .par_iter()
                    .map(|tile| {
                        let filter_path = run_path.join(format!(
                            "Data/Intensities/BaseCalls/L{:03}/s_{}_{}.filter",
                            lane, lane, tile,
                        ));
                        let filter = match filter_decoder(&filter_path) {
                            Ok(filter) => filter,
                            Err(e) => {
                                panic!("Error reading filter {} {}", filter_path.display(), e)
                            }
                        };

                        (tile.clone(), filter)
                    })
                    .unzip_into_vecs(&mut lane_surface_tile_ids, &mut lane_surface_filters);

                let mut lane_surface_n_pfs = Vec::new();
                let mut lane_surface_pf_filters = Vec::new();

                lane_surface_filters
                    .par_iter()
                    .map(|filter| {
                        let n_pf: usize = filter.iter().map(|&b| [0, 1, 1, 2][b as usize]).sum();

                        let mut pf_filter: Vec<_> = std::iter::repeat(3).take(n_pf / 2).collect();
                        if n_pf % 2 == 1 {
                            pf_filter.push(2)
                        }

                        (n_pf, pf_filter)
                    })
                    .unzip_into_vecs(&mut lane_surface_n_pfs, &mut lane_surface_pf_filters);

                info!("loaded {} filters and ids", lane_surface_filters.len());

                filters.insert([lane, surface], lane_surface_filters);
                pf_filters.insert([lane, surface], lane_surface_pf_filters);
                tile_ids.insert([lane, surface], lane_surface_tile_ids);
                n_pfs.insert([lane, surface], lane_surface_n_pfs);

                read_headers.insert([lane, surface], lane_surface_read_headers);
                index_headers.insert([lane, surface], lane_surface_index_headers);
            }
        }

        // check to make sure our "constant qscore map" assumption is correct
        let mut qscore_maps: std::collections::HashSet<_> = index_headers
            .values()
            .flat_map(|hs| hs.iter().flat_map(|h| h.iter().map(|h| h.bins.clone())))
            .collect();

        qscore_maps.extend(
            read_headers
                .values()
                .flat_map(|hs| hs.iter().flat_map(|h| h.iter().map(|h| h.bins.clone()))),
        );

        if qscore_maps.len() == 1 {
            if !qscore_maps.contains(&vec![35, 44, 58, 70]) {
                panic!("Got a new qscore map! {:?}", qscore_maps);
            }
        } else {
            panic!("Got {} different qscore maps!", qscore_maps.len());
        }

        debug!("qscore map: {:?}", qscore_maps);

        let novaseq_run = NovaSeqRun {
            run_path,
            run_info,
            run_id,
            locs,
            filters,
            pf_filters,
            tile_ids,
            n_pfs,
            read_headers,
            index_headers,
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
        let novaseq_run = NovaSeqRun::read_path(run_path, false).unwrap();

        assert_eq!(novaseq_run.run_id, "@A00111:296:HJCWWDSXX");
    }

    #[test]
    fn test_index_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path, true).unwrap();

        assert_eq!(novaseq_run.run_id, "@A00111:296:HJCWWDSXX");
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn no_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXXX");
        NovaSeqRun::read_path(run_path, false).unwrap();
    }
}
