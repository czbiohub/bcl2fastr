//! Represents a NovaSeq sequencing run as a struct 
//! that can be shared across threads

use std::path::PathBuf;
use std::collections::HashMap;

use ndarray::Array2;
use rayon::prelude::*;

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
    /// a string with the run info formatted for read headers
    pub run_id: String,
    /// start and end indices for the reads
    pub read_ix: Vec<(usize, usize)>,
    /// start and end indices for the indexes
    pub index_ix: Vec<(usize, usize)>,
    /// a single universal locs array, same for every tile
    pub locs: Locs,
    /// a map from (lane, surface) tuples to vectors of filters
    pub filters: HashMap<(u32, u32), Vec<Filter>>,
    /// a map from (lane, surface) tuples to vectors of post-filter "filters"
    /// which are needed to remove empty half-bytes at the end of tiles
    pub pf_filters: HashMap<(u32, u32), Vec<Filter>>,
    /// a map from (lane, surface) tuples to vectors of indexes into locs
    /// for the combined filters
    pub read_ids: HashMap<(u32, u32), Vec<Array2<u32>>>,
    /// a map from (lane, surface) tuples to vectors of CBCL headers
    pub headers: HashMap<(u32, u32), Vec<CBCLHeader>>,
}


impl NovaSeqRun {
    /// Loads a NovaSeqRun located in `run_path`, aggregating filters and headers
    /// in chunks of `tile_chunk` tiles each. If `index_only` is true, will only
    /// load in data for cycles that are in indexes, and adjusts `index_ix` attribute
    /// accordingly. Uses threads to load the data in parallel.
    pub fn read_path(
        run_path: PathBuf, tile_chunk: usize, index_only: bool,
    ) -> std::io::Result<NovaSeqRun> {
        let run_info = parse_run_info(&run_path.join("RunInfo.xml"))?;
        let run_id = format!(
            "@{}:{}:{}",
            run_info.instrument,
            run_info.number,
            run_info.flowcell,
        );

        let locs = locs_decoder(&run_path.join("Data/Intensities/s.locs"))?;

        let mut headers = HashMap::new();
        let mut filters = HashMap::new();
        let mut pf_filters = HashMap::new();
        let mut read_ids = HashMap::new();

        let n_lanes = run_info.flowcell_layout.lane_count;
        let n_surfaces = run_info.flowcell_layout.surface_count;
        let n_cycles = run_info.reads.iter().map(|v| v.num_cycles).sum();

        let read_ix: Vec<(usize, usize)> = if index_only { 
            Vec::new() 
        } else {
            run_info.reads.iter()
                .filter(|r| !r.is_indexed_read)
                .map(|r| (r.start, r.end))
                .collect()
        };

        let indexes: Vec<_> = run_info.reads.iter()
            .filter(|r| r.is_indexed_read)
            .collect();

        let index_ix: Vec<(usize, usize)> = if index_only {
            indexes.iter().scan(0, |i, r| {
                *i += r.num_cycles;
                Some((*i - r.num_cycles, *i))
            }).collect()
        } else {
            indexes.iter().map(|r| (r.start, r.end)).collect()
        };

        let cycle_range = if index_only {
            (indexes[0].start + 1)..=(indexes[indexes.len() - 1].end)
        } else {
            1..=n_cycles
        };

        for lane in 1..=n_lanes {
            for surface in 1..=n_surfaces {
                println!("lane {} - surface {}", lane, surface);
                let lane_surface_headers: Vec<CBCLHeader> = cycle_range.clone()
                    .into_par_iter()
                    .map( |cycle| {
                        let cbcl_path = run_path.join(
                            format!(
                                "Data/Intensities/BaseCalls/L{:03}/C{}.1/L{:03}_{}.cbcl",
                                lane,
                                cycle,
                                lane,
                                surface,
                            )
                        );

                        cbcl_header_decoder(&cbcl_path, tile_chunk).unwrap()
                    }
                ).collect();

                println!("loaded {} headers", lane_surface_headers.len());

                // tile numbers are not stored by surface in RunInfo, so we are
                // taking advantage of the headers having the right names
                let filter_and_id_vec: Vec<(_, _, _)> = lane_surface_headers[0].tiles
                    .par_chunks(tile_chunk)
                    .map( |tile_chunk| {
                        // when the tiles are not filtered, we need a combined filter 
                        // for all of the tiles being extracted
                        let mut filter = Filter::new();
                        // when the tiles are already filtered, we need to take into 
                        // account any half-packed bytes at the end of each filter
                        let mut pf_filter = Filter::new();
                        let mut read_id: Vec<[u32; 3]> = Vec::new();

                        for tile in tile_chunk {
                            let filter_path = run_path.join(
                                format!(
                                    "Data/Intensities/BaseCalls/L{:03}/s_{}_{}.filter",
                                    lane,
                                    lane,
                                    tile,
                                )
                            );
                            let tile_filter = filter_decoder(&filter_path).unwrap();

                            read_id.extend(tile_filter.iter().enumerate().filter_map(
                                |(i, &b)|
                                if b {
                                    Some([tile.clone(), locs[i][0], locs[i][1]])
                                } else {
                                    None
                                }
                            ));

                            let n_pf = tile_filter.iter().map(
                                |&b| if b { 1 } else { 0 }
                            ).sum::<usize>();

                            pf_filter.extend(std::iter::repeat(true).take(n_pf));
                            if n_pf % 2 == 1 {
                                pf_filter.push(false)
                            }

                            filter.extend(tile_filter);
                        }

                        // store read_id information consisting of tile number & xy locs
                        let read_id: Array2<u32> = match Array2::from_shape_vec(
                            (read_id.len(), 3), read_id.iter().flatten().cloned().collect()
                        ) {
                            Ok(a) => a,
                            Err(e) => panic!("Failed to create read_id array: {}", e)
                        };

                        (filter, pf_filter, read_id)
                    })
                    .collect();

                let mut lane_surface_filters = Vec::new();
                let mut lane_surface_pf_filters = Vec::new();
                let mut lane_surface_read_ids = Vec::new();

                for (filter, pf_filter, read_id) in filter_and_id_vec {
                    lane_surface_filters.push(filter);
                    lane_surface_pf_filters.push(pf_filter);
                    lane_surface_read_ids.push(read_id);
                }

                println!("loaded {} filters and ids", lane_surface_filters.len());

                headers.insert((lane, surface), lane_surface_headers);
                filters.insert((lane, surface), lane_surface_filters);
                pf_filters.insert((lane, surface), lane_surface_pf_filters);
                read_ids.insert((lane, surface), lane_surface_read_ids);
            }
        }

        let novaseq_run = NovaSeqRun {
            run_path,
            run_info,
            run_id,
            read_ix,
            index_ix,
            locs,
            filters,
            pf_filters,
            read_ids,
            headers,
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
        let expected_read_ix = vec![(0, 4), (20, 24)];
        let expected_index_ix = vec![(4, 12), (12, 20)];
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();

        assert_eq!(novaseq_run.run_id, "@A00111:296:HJCWWDSXX");
        assert_eq!(novaseq_run.read_ix, expected_read_ix);
        assert_eq!(novaseq_run.index_ix, expected_index_ix);
    }

    #[test]
    fn test_index_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let expected_read_ix = vec![];
        let expected_index_ix = vec![(0, 8), (8, 16)];
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, true).unwrap();

        assert_eq!(novaseq_run.run_id, "@A00111:296:HJCWWDSXX");
        assert_eq!(novaseq_run.read_ix, expected_read_ix);
        assert_eq!(novaseq_run.index_ix, expected_index_ix);
    }

    #[test]
    #[should_panic(
        expected = r#"No such file or directory"#
    )]
    fn no_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXXX");
        NovaSeqRun::read_path(run_path, 2, false).unwrap();
    }
}
