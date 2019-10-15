//! Extract and decompress a set of tiles from a vector of cbcl files.

extern crate flate2;
extern crate ndarray;

use std::{
    fs::File,
    io::prelude::*,
    io::SeekFrom,
};
use flate2::read::MultiGzDecoder;

use crate::cbcl_header_decoder::CBCLHeader;
use crate::filter_decoder::Filter;
use crate::novaseq_run::NovaSeqRun;


/// unpacks a single byte into four 2-bit integers
fn unpack_byte(b: &u8) -> Vec<(u8, u8)> {
    let q_1 = (b >> 6) & 3u8;
    let b_1 = (b >> 4) & 3u8;
    let q_2 = (b >> 2) & 3u8;
    let b_2 = b & 3u8;

    vec![(b_2, q_2), (b_1, q_1)]
}


/// extract multiple tiles from a CBCL file and return decompressed bytes
fn extract_tiles(
    header: &CBCLHeader, first_idx: usize, last_idx: usize
) -> std::io::Result<Vec<u8>> {
    // calculate start byte position
    let start_pos = (
        header.header_size 
        + header.tile_offsets[0..first_idx]
                .iter()
                .map(|v| v[3])
                .sum::<u32>()
    ) as u64;

    // calculate end byte position and expected size of the decompressed tile(s)
    // index 2 is the uncompressed tile size
    let uncompressed_size = header.tile_offsets[first_idx..last_idx]
                                  .iter()
                                  .map(|v| v[2])
                                  .sum::<u32>() as usize;
    // index 3 is the compressed tile size
    let compressed_size = header.tile_offsets[first_idx..last_idx]
                                .iter()
                                .map(|v| v[3])
                                .sum::<u32>() as usize;

    // open file and read whole file into a buffer
    let mut cbcl = File::open(&header.cbcl_path)?;
    cbcl.seek(SeekFrom::Start(start_pos))?;

    let mut read_buffer = vec![0u8; compressed_size];
    cbcl.read_exact(&mut read_buffer)?;

    // use MultiGzDecoder to uncompress the number of bytes summed 
    // over the offsets of all tile_idces
    let mut uncomp_bytes = vec![0u8; uncompressed_size];
    let mut gz = MultiGzDecoder::new(&read_buffer[..]);
    gz.read_exact(&mut uncomp_bytes)?;

    Ok(uncomp_bytes)
}


/// given a CBCL file and some tiles: extract, translate and filter the bases+scores
fn process_tiles(
    header: &CBCLHeader, filter: &Filter, first_idx: usize, last_idx: usize
) -> std::io::Result<Vec<(u8, u8)>> {
    let uncomp_bytes = extract_tiles(header, first_idx, last_idx)?;

    // unpack the bytes into tuples (two per byte), then use the filter to filter 
    let bq_pairs = uncomp_bytes.iter()
                               .map(|v| unpack_byte(v))
                               .flatten()
                               .zip(filter)
                               .filter_map(|v| if *v.1 { Some(v.0) } else { None })
                               .collect();

    Ok(bq_pairs)
}


/// Create arrays of read and qscore values from a set of tiles
pub fn extract_reads(
    headers: &Vec<CBCLHeader>, filters: &Vec<Filter>, first_idx: usize, last_idx: usize
) -> std::io::Result<(ndarray::Array2<u8>, ndarray::Array2<u8>)> {
    // TODO: calculate filters in startup

    // when the tiles are not filtered, we need a combined filter for
    // all of the tiles being extracted
    let combined_filter: Vec<bool> = filters[first_idx..last_idx].iter()
                                                                 .cloned()
                                                                 .flatten()
                                                                 .collect();

    // when the tiles are already filtered, we need to take into account
    // any half-packed bytes at the end
    let mut ex_filter = Vec::new();

    for filter in filters[first_idx..last_idx].iter() {
        let n_pf = filter.iter().map(|&b| if b { 1 } else { 0 }).sum::<u64>();
        ex_filter.append(&mut vec![true; n_pf as usize]);
        if n_pf % 2 == 1 {
            ex_filter.push(false)
        }
    }

    let n_pf = ex_filter.iter().map(|&b| if b { 1 } else { 0 }).sum::<u64>() as usize;
    let n_cycles = headers.len();
    let mut read_array = ndarray::Array2::from_elem((n_cycles, n_pf), 4u8);
    let mut qscore_array = ndarray::Array2::from_elem((n_cycles, n_pf), 4u8);

    for (i, h) in headers.iter().enumerate() {
        let mut read_row = read_array.index_axis_mut(ndarray::Axis(0), i);
        let mut qscore_row = qscore_array.index_axis_mut(ndarray::Axis(0), i);

        if h.non_pf_clusters_excluded {
            if let Ok(tile_bytes) = process_tiles(
                &h, &ex_filter, first_idx, last_idx
            ) {
                let (b_array, q_array): (Vec<u8>, Vec<u8>) = tile_bytes.iter()
                                                                       .cloned()
                                                                       .unzip();
                read_row.assign(&ndarray::ArrayView::from(&b_array));
                qscore_row.assign(&ndarray::ArrayView::from(&q_array));
            };
        } else {
            if let Ok(tile_bytes) = process_tiles(
                &h, &combined_filter, first_idx, last_idx
            ) {
                let (b_array, q_array): (Vec<u8>, Vec<u8>) = tile_bytes.iter()
                                                                       .cloned()
                                                                       .unzip();
                read_row.assign(&ndarray::ArrayView::from(&b_array));
                qscore_row.assign(&ndarray::ArrayView::from(&q_array));
            };
        }
    }

    Ok((read_array, qscore_array))
}


/// Iterate through all lanes and surfaces of a run and extract tiles in chunks
pub fn extract_samples(
    novaseq_run: NovaSeqRun, tile_chunk: usize
) -> Result<(), &'static str> {
    for lane in 1..=novaseq_run.run_info.runs.flow_cell_layout.lane_count {
        for surface in 1..=novaseq_run.run_info.runs.flow_cell_layout.surface_count {
            let headers = novaseq_run.headers.get(&(lane, surface)).unwrap();
            let filters = novaseq_run.filters.get(&(lane, surface)).unwrap();

            let n_tiles = headers[0].tile_offsets.len();

            for first_idx in (0..n_tiles).step_by(tile_chunk) {
                let last_idx = std::cmp::min(first_idx + tile_chunk, n_tiles);

                let res = extract_reads(&headers, &filters, first_idx, last_idx);
                println!("extracted some reads! {:?}", res);
            }
        }
    }

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use crate::cbcl_header_decoder::cbcl_header_decoder;

    #[test]
    fn extract_tiles() {
        let cbcl_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX").join(
            "Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl"
        );
        let cbcl_header = cbcl_header_decoder(&cbcl_path).unwrap();

        let expected_bytes = vec![
            212, 254, 220, 221, 166, 108, 217, 232, 236, 221,
            157, 216, 220, 220, 205, 222, 140, 212, 157, 254,
            199, 221, 237, 185, 252, 199, 237, 253, 253, 68,
            237, 205, 199, 199, 237, 109, 205, 79, 200, 220,
            76, 253, 204, 253, 95, 223, 238, 78, 79, 206,
            220, 152, 220, 157, 255, 196, 207, 207, 133, 78,
            236, 222, 205, 254, 237, 204, 198, 218, 236, 204,
            206, 204, 214, 207, 222, 204, 201, 221, 103, 207,
            204, 196, 204, 88, 216, 205, 222, 251, 253, 206,
            206, 237, 223, 220, 205, 76, 220, 205, 232, 220
        ];

        let uncomp_bytes = super::extract_tiles(&cbcl_header, 0, 2).unwrap();

        assert_eq!(uncomp_bytes, expected_bytes)
    }

    #[test]
    fn process_tiles() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path).unwrap();
        let first_idx = 0;
        let last_idx = 1;

        let expected_bq_pairs = vec![
            (3, 3), (0, 3), (1, 3), (1, 3), (1, 3), (0, 3), (1, 2), (1, 3)
        ];

        let header = &novaseq_run.headers.get(&(1, 1)).unwrap()[0];
        let combined_filter = novaseq_run.filters.get(&(1, 1))
                                                 .unwrap()[first_idx..last_idx]
                                                 .iter()
                                                 .cloned()
                                                 .flatten()
                                                 .collect();

        let bq_pairs: Vec<(u8, u8)> = super::process_tiles(
            header,
            &combined_filter,
            0, 1,
        ).unwrap().into_iter().take(8).collect();

        assert_eq!(bq_pairs, expected_bq_pairs)
    }

    #[test]
    fn extract_samples() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path).unwrap();

        super::extract_samples(novaseq_run, 2).unwrap()
    }
}
