//! Represents a NovaSeq sequencing run as a struct
//! that can be shared across threads

use std::collections::HashMap;
use std::path::PathBuf;

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
    /// pre-compute this size once: the maximum number of bytes in any block,
    /// used to preallocate a vector while reading the data
    pub max_vec_size: usize,
    /// a single universal locs array, same for every tile
    pub locs: Locs,
    /// a map from [lane, surface] tuples to vectors of filters
    pub filters: HashMap<[usize; 2], Vec<Filter>>,
    /// a map from [lane, surface] tuples to vectors of post-filter "filters"
    /// which are needed to remove empty half-bytes at the end of tiles
    pub pf_filters: HashMap<[usize; 2], Vec<Filter>>,
    /// a map from [lane, surface] to vectors of tiles, so we can produce
    /// the read header
    pub tile_ids: HashMap<[usize; 2], Vec<Vec<(u32, usize)>>>,
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
    pub fn read_path(
        run_path: PathBuf,
        tile_chunk: usize,
        index_only: bool,
    ) -> std::io::Result<NovaSeqRun> {
        let run_info = parse_run_info(&run_path.join("RunInfo.xml"))?;
        let run_id = format!(
            "@{}:{}:A{}",
            run_info.instrument, run_info.number, run_info.flowcell,
        );

        // need to repeat locs for each tile in tile_chunk
        let locs = locs_decoder(&run_path.join("Data/Intensities/s.locs"))?;
        let locs = locs
            .iter()
            .cycle()
            .take(locs.len() * tile_chunk)
            .cloned()
            .collect();

        let mut read_headers = HashMap::new();
        let mut index_headers = HashMap::new();

        let mut filters = HashMap::new();
        let mut pf_filters = HashMap::new();
        let mut tile_ids = HashMap::new();

        let mut max_vec_size: Vec<usize> = Vec::new();

        for lane in 1..=run_info.flowcell_layout.lane_count {
            for surface in 1..=run_info.flowcell_layout.surface_count {
                println!("lane {} - surface {}", lane, surface);

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
                                format!("L{:03}/C{}.1/L{:03}_{}.cbcl", lane, cycle, lane, surface,),
                            );

                            CBCLHeader::from_path(&cbcl_path, tile_chunk).unwrap()
                        })
                        .collect();

                    max_vec_size.push(
                        these_headers
                            .iter()
                            .flat_map(|h| h.uncompressed_size.iter().cloned())
                            .max()
                            .unwrap(),
                    );

                    if read.is_indexed_read {
                        lane_surface_index_headers.push(these_headers);
                    } else {
                        lane_surface_read_headers.push(these_headers);
                    }
                }

                // will always have index_headers, not always headers

                // tile numbers are not stored by surface in RunInfo, so we are
                // taking advantage of the headers having the right names
                let filter_and_id_vec = lane_surface_index_headers[0][0]
                    .tiles
                    .par_chunks(tile_chunk)
                    .map(|tile_chunk| {
                        // when the tiles are not filtered, we need a combined filter
                        // for all of the tiles being extracted
                        let mut filter = Filter::new();
                        // when the tiles are already filtered, we need to take into
                        // account any half-packed bytes at the end of each filter
                        let mut pf_filter = Filter::new();
                        let mut tile_id = Vec::new();

                        for tile in tile_chunk {
                            let filter_path = run_path.join(format!(
                                "Data/Intensities/BaseCalls/L{:03}/s_{}_{}.filter",
                                lane, lane, tile,
                            ));
                            let tile_filter = filter_decoder(&filter_path).unwrap();

                            let n_pf = tile_filter.iter().map(|&b| if b { 1 } else { 0 }).sum();

                            tile_id.push((tile.clone(), n_pf));

                            pf_filter.extend(std::iter::repeat(true).take(n_pf));
                            if n_pf % 2 == 1 {
                                pf_filter.push(false)
                            }

                            filter.extend(tile_filter);
                        }

                        (filter, pf_filter, tile_id)
                    })
                    .collect::<Vec<_>>();

                let mut lane_surface_filters = Vec::new();
                let mut lane_surface_pf_filters = Vec::new();
                let mut lane_surface_tile_ids = Vec::new();

                for (filter, pf_filter, tile_id) in filter_and_id_vec {
                    lane_surface_filters.push(filter);
                    lane_surface_pf_filters.push(pf_filter);
                    lane_surface_tile_ids.push(tile_id);
                }

                println!("loaded {} filters and ids", lane_surface_filters.len());

                filters.insert([lane, surface], lane_surface_filters);
                pf_filters.insert([lane, surface], lane_surface_pf_filters);
                tile_ids.insert([lane, surface], lane_surface_tile_ids);
                read_headers.insert([lane, surface], lane_surface_read_headers);
                index_headers.insert([lane, surface], lane_surface_index_headers);
            }
        }

        let max_vec_size = max_vec_size.iter().cloned().max().unwrap();
        // round up to nearest multiple of 8
        let max_vec_size = if max_vec_size % 8 > 0 {
            max_vec_size + 8 - (max_vec_size % 8)
        } else {
            max_vec_size
        };

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

        let novaseq_run = NovaSeqRun {
            run_path,
            run_info,
            run_id,
            max_vec_size,
            locs,
            filters,
            pf_filters,
            tile_ids,
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
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();

        assert_eq!(novaseq_run.run_id, "@A00111:296:AHJCWWDSXX");
    }

    #[test]
    fn test_index_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, true).unwrap();

        assert_eq!(novaseq_run.run_id, "@A00111:296:AHJCWWDSXX");
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn no_run() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXXX");
        NovaSeqRun::read_path(run_path, 2, false).unwrap();
    }
}
