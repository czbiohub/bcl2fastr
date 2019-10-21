//! Extract and decompress a set of tiles from a vector of cbcl files.

use std::{
    fs::File,
    io::prelude::*,
    io::SeekFrom,
    path::PathBuf
};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use ndarray::{Array2, ArrayView, Axis};

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


/// converts from 0..3 values to the appropriate base, or N if the qscore is too low
fn u8_to_base(b: u8, q: u8) -> u8 {
    if q <= 35 {
        return b'N'
    }

    match b {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N'
    }
}


/// extract multiple tiles from a CBCL file and return decompressed bytes
fn extract_tiles(header: &CBCLHeader, i: usize) -> std::io::Result<Vec<u8>> {
    let start_pos = header.start_pos[i];
    let uncompressed_size = header.uncompressed_size[i];
    let compressed_size = header.compressed_size[i];

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
    header: &CBCLHeader, filter: &Filter, i: usize,
) -> std::io::Result<Vec<(u8, u8)>> {
    let uncomp_bytes = extract_tiles(header, i)?;

    // unpack the bytes into tuples (two per byte), then use the filter to filter 
    let bq_pairs = uncomp_bytes.iter()
                               .map(|v| unpack_byte(v))
                               .flatten()
                               .zip(filter)
                               .filter_map(|v| if v.1 { Some(v.0) } else { None })
                               .collect();

    Ok(bq_pairs)
}


/// Create arrays of read and qscore values from a set of tiles
pub fn extract_reads(
    headers: &Vec<CBCLHeader>, filter: &Filter, pf_filter: &Filter, i: usize,
) -> std::io::Result<(Array2<u8>, Array2<u8>)> {
    let n_pf = pf_filter.iter().map(|b| if b { 1 } else { 0 }).sum::<u64>() as usize;
    let n_cycles = headers.len();

    let mut read_array = Array2::from_elem((n_cycles, n_pf), 4u8);
    let mut qscore_array = Array2::from_elem((n_cycles, n_pf), 4u8);

    for (j, h) in headers.iter().enumerate() {
        let mut read_row = read_array.index_axis_mut(Axis(0), j);
        let mut qscore_row = qscore_array.index_axis_mut(Axis(0), j);

        let h_filter = if h.non_pf_clusters_excluded { pf_filter } else { filter };

        if let Ok(tile_bytes) = process_tiles(&h, &h_filter, i) {
            let (b_array, q_array): (Vec<u8>, Vec<u8>) = 
                tile_bytes.iter()
                          .cloned()
                          .map(|(b, q)| (b, h.decode_qscore(q)))
                          .map(|(b, q)| (u8_to_base(b, q), q))
                          .unzip();

            read_row.assign(&ArrayView::from(&b_array));
            qscore_row.assign(&ArrayView::from(&q_array));
        };
    }

    Ok((read_array, qscore_array))
}


/// Iterate through all lanes and surfaces of a run and extract tiles in chunks
pub fn extract_samples(
    novaseq_run: NovaSeqRun, output_path: PathBuf
) -> Result<(), &'static str> {
    let out_file = match File::create(output_path.join("Undetermined.fastq.gz")) {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    let mut gz_writer = GzEncoder::new(out_file, flate2::Compression::new(9));

    for lane in 1..=novaseq_run.run_info.runs.flow_cell_layout.lane_count {
        for surface in 1..=novaseq_run.run_info.runs.flow_cell_layout.surface_count {
            let headers = novaseq_run.headers.get(&(lane, surface)).unwrap();
            let filters = novaseq_run.filters.get(&(lane, surface)).unwrap();
            let pf_filters = novaseq_run.pf_filters.get(&(lane, surface)).unwrap();
            let read_ids = novaseq_run.read_ids.get(&(lane, surface)).unwrap();

            for (i, ((filter, pf_filter), rid)) in filters.iter()
                                                          .zip(pf_filters)
                                                          .zip(read_ids)
                                                          .enumerate() {
                if let Ok((reads, qscores)) = extract_reads(&headers, filter, pf_filter, i) {
                    for ((r_row, q_row), id_row) in reads.axis_iter(Axis(1))
                                                         .zip(qscores.axis_iter(Axis(1)))
                                                         .zip(rid.axis_iter(Axis(0))) {

                        write!(
                            &mut gz_writer,
                            "{}:{}:{}:{}:{}\n",
                            novaseq_run.run_id,
                            lane,
                            id_row[0],
                            id_row[1],
                            id_row[2],
                        ).unwrap();

                        gz_writer.write(&r_row.to_vec()[..]).unwrap();
                        gz_writer.write(b"\n+\n").unwrap();
                        gz_writer.write(&q_row.to_vec()[..]).unwrap();
                        gz_writer.write(b"\n").unwrap();
                    }
                };
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
        let cbcl_header = cbcl_header_decoder(&cbcl_path, 2).unwrap();

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

        let uncomp_bytes = super::extract_tiles(&cbcl_header, 0).unwrap();

        assert_eq!(uncomp_bytes, expected_bytes)
    }

    #[test]
    fn process_tiles() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path, 2).unwrap();

        let expected_bq_pairs = vec![
            (3, 3), (0, 3), (1, 3), (1, 3), (1, 3), (0, 3), (1, 2), (1, 3)
        ];

        let header = &novaseq_run.headers.get(&(1, 1)).unwrap()[0];
        let filter = &novaseq_run.filters.get(&(1, 1)).unwrap()[0];

        let bq_pairs: Vec<(u8, u8)> = super::process_tiles(
            header, filter, 0
        ).unwrap().into_iter().take(8).collect();

        assert_eq!(bq_pairs, expected_bq_pairs)
    }
}
