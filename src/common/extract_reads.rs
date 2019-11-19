//! Extract and decompress a set of tiles from a vector of cbcl files.

use std::{
    fs::File,
    io::prelude::*,
    io::SeekFrom,
};

use flate2::read::MultiGzDecoder;
use ndarray::{Array3, ArrayViewMut2, Axis, ShapeBuilder};

use crate::cbcl_header_decoder::CBCLHeader;


/// converts from 0..3 values to the appropriate base, or N if the qscore is too low
fn u8_to_base(b: u8, q: u8) -> u8 {
    if q <= 35 { return b'N' }

    match b {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}


/// unpacks a single byte into four 2-bit integers
fn unpack_byte(b: &u8, filter: &[bool], header: &CBCLHeader) -> Vec<u8> {
    match filter {
        [true, true] => {
            let q_1 = header.decode_qscore((b >> 6) & 3u8);
            let b_1 = u8_to_base((b >> 4) & 3u8, q_1);
            let q_2 = header.decode_qscore((b >> 2) & 3u8);
            let b_2 = u8_to_base(b & 3u8, q_2);
        
            vec![b_2, q_2, b_1, q_1]
        },
        [true, false] => {
            let q_2 = header.decode_qscore((b >> 2) & 3u8);
            let b_2 = u8_to_base(b & 3u8, q_2);
        
            vec![b_2, q_2]
        },
        [false, true] => {
            let q_1 = header.decode_qscore((b >> 6) & 3u8);
            let b_1 = u8_to_base((b >> 4) & 3u8, q_1);
            vec![b_1, q_1]
        },
        _ => vec![],
    }
}


/// extract multiple tiles from a CBCL file and return decompressed bytes
fn extract_tiles(
    header: &CBCLHeader, i: usize, read_buffer: &mut Vec<u8>
) -> std::io::Result<()> {
    let start_pos = header.start_pos[i];
    let uncompressed_size = header.uncompressed_size[i];

    // open file and seek to start position
    let mut cbcl = File::open(&header.cbcl_path)?;
    cbcl.seek(SeekFrom::Start(start_pos))?;

    // use MultiGzDecoder to decompress
    let mut gz = MultiGzDecoder::new(&cbcl);
    gz.read_exact(&mut read_buffer[..uncompressed_size])?;

    Ok(())
}


/// given a CBCL file and some tiles: extract, translate and filter the bases+scores
fn process_tiles(
    read_buffer: &mut Vec<u8>,
    bq_cycle: &mut ArrayViewMut2<u8>,
    header: &CBCLHeader,
    filter: &[bool],
    i: usize,
) -> () {
    if let Ok(_) = extract_tiles(header, i, read_buffer) {
        // unpack the bytes, filtering out the reads that didn't pass
        bq_cycle.iter_mut().zip(
            read_buffer.iter()
                .zip(filter.chunks(2))
                .flat_map(|(v, f)| unpack_byte(v, f, header))
        ).for_each(|(a, b)| { *a = b; });
    } else {
        bq_cycle.index_axis_mut(Axis(1), 0).fill(b'N');
        bq_cycle.index_axis_mut(Axis(1), 1).fill(b'#');
    }
}


/// Create arrays of read and qscore values from a set of tiles
pub fn extract_read_block(
    headers: &[CBCLHeader],
    filter: &[bool],
    pf_filter: &[bool],
    rq_array: &mut Array3<u8>,
    read_buffer: &mut Vec<u8>,
    i: usize,
    n_pf: usize,
) -> () {
    let n_cycles = headers.len();
    let mut bq_array = rq_array.slice_mut(ndarray::s![..n_cycles, ..n_pf, ..]);

    for (mut cycle, h) in bq_array.axis_iter_mut(Axis(0)).zip(headers) {
        let h_filter = if h.non_pf_clusters_excluded { pf_filter } else { filter };

        process_tiles(read_buffer, &mut cycle, h, h_filter, i);
    }
}


/// given a vector of index headers, extract the indices as byte vectors, then re-zip
/// them into the groups of indices corresponding to each read
pub fn extract_indices<'a>(
    headers: &[Vec<CBCLHeader>],
    filter: &[bool],
    pf_filter: &[bool],
    read_buffer: &mut Vec<u8>,
    i: usize,
    n_pf: usize,
) -> Vec<Array3<u8>> {
    let mut index_arrays: Vec<_> = headers.iter().map( |h| 
        Array3::zeros((h.len(), n_pf, 2).f())
    ).collect();

    headers.iter()
        .zip(index_arrays.iter_mut())
        .for_each( |(h, mut out_arr)|
            extract_read_block(h, filter, pf_filter, &mut out_arr, read_buffer, i, n_pf)
        );

    index_arrays
}


#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use ndarray::Array2;

    use crate::cbcl_header_decoder::cbcl_header_decoder;
    use crate::novaseq_run::NovaSeqRun;

    #[test]
    fn u8_to_base() {
        let expected_bases = vec![b'A', b'C', b'G', b'T', b'N'];
        let actual_bases: Vec<_> = [0, 1, 2, 3, 4].iter()
            .map(|&b| super::u8_to_base(b, 70))
            .collect();

        assert_eq!(actual_bases, expected_bases);
    }

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

        let mut read_buffer = vec![0u8; cbcl_header.uncompressed_size[0]];

        super::extract_tiles(&cbcl_header, 0, &mut read_buffer).unwrap();

        assert_eq!(read_buffer, expected_bytes)
    }

    #[test]
    fn process_tiles() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();

        let expected_bq_pairs = vec![
            84, 70, 65, 70, 67, 70, 67, 70, 67, 70, 65, 70, 67, 58, 67, 70
        ];

        let header = &novaseq_run.read_headers.get(&[1, 1]).unwrap()[0][0];
        let filter = &novaseq_run.filters.get(&[1, 1]).unwrap()[0];

        let n_pf = filter.iter().map(|&b| if b { 1 } else { 0 }).sum();
        let mut bq_array = Array2::zeros((n_pf, 2));
        let mut read_buffer = vec![0u8; novaseq_run.max_vec_size];

        super::process_tiles(
            &mut read_buffer, &mut bq_array.view_mut(), header, filter, 0
        );

        let bq_pairs: Vec<_> = bq_array.iter().cloned().take(16).collect();

        assert_eq!(bq_pairs, expected_bq_pairs)
    }
}
