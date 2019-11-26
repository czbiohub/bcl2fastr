//! Extract and decompress a set of tiles from a vector of cbcl files.

use std::{fs::File, io::prelude::*, io::SeekFrom};

use bitpacking::{BitPacker, BitPacker1x};
use flate2::read::MultiGzDecoder;
use ndarray::{ArrayViewMut3, Axis};

use crate::cbcl_header_decoder::CBCLHeader;

/// converts from 0..3 values to a PHRED qscore. Assumes a constant mapping
#[inline]
fn u8_to_qscore(b: u8) -> u8 {
    match b {
        0 => 35,
        1 => 44,
        2 => 58,
        3 => 70,
        _ => b'#',
    }
}

/// converts from 0..3 values to the appropriate base, or N if the qscore is too low
#[inline]
fn u8_to_base(b: u8, q: u8) -> u8 {
    if q <= 35 {
        return b'N';
    } else {
        match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        }
    }
}

/// extract multiple tiles from a CBCL file and return decompressed bytes
fn extract_tiles(
    header: &CBCLHeader,
    tile_i: usize,
    read_buffer: &mut Vec<u8>,
) -> std::io::Result<()> {
    let start_pos = header.start_pos[tile_i];
    let uncompressed_size = header.uncompressed_size[tile_i];

    // open file and seek to start position
    let mut cbcl = File::open(&header.cbcl_path)?;
    cbcl.seek(SeekFrom::Start(start_pos))?;

    // use MultiGzDecoder to decompress
    let mut gz = MultiGzDecoder::new(&cbcl);
    gz.read_exact(&mut read_buffer[..uncompressed_size])?;

    Ok(())
}

/// Extract, translate and filter the bases+scores into a pre-allocated array
pub fn extract_read_block(
    headers: &[CBCLHeader],
    filter: &[bool],
    pf_filter: &[bool],
    bq_array: &mut ArrayViewMut3<u8>,
    read_buffer: &mut Vec<u8>,
    decompression_buffer: &mut Vec<u32>,
    tile_i: usize,
) -> () {
    let bitpacker = BitPacker1x::new();

    for (mut bq_cycle, header) in bq_array.axis_iter_mut(Axis(0)).zip(headers) {
        let h_filter = if header.non_pf_clusters_excluded {
            pf_filter
        } else {
            filter
        };

        if let Ok(_) = extract_tiles(header, tile_i, read_buffer) {
            // unpack the bytes
            read_buffer
                .chunks(8)
                .zip(decompression_buffer.chunks_exact_mut(32))
                .for_each(|(b, mut d)| {
                    bitpacker.decompress(b, &mut d, 2);
                });

            bq_cycle
                .iter_mut()
                .zip(
                    decompression_buffer
                        .chunks(2)
                        .zip(h_filter)
                        .filter_map(|(v, &b)| if b { Some(v) } else { None })
                        .flatten(),
                )
                .for_each(|(a, b)| {
                    *a = *b as u8;
                });

            bq_cycle.axis_iter_mut(Axis(0)).for_each(|mut rq| {
                rq[1] = u8_to_qscore(rq[1]);
                rq[0] = u8_to_base(rq[0], rq[1]);
            });
        } else {
            bq_cycle.index_axis_mut(Axis(1), 0).fill(b'N');
            bq_cycle.index_axis_mut(Axis(1), 1).fill(b'#');
        }
    }
}

#[cfg(test)]
mod tests {
    use ndarray::{Array3, ShapeBuilder};
    use std::path::PathBuf;

    use crate::cbcl_header_decoder::CBCLHeader;
    use crate::novaseq_run::NovaSeqRun;

    #[test]
    fn u8_to_base() {
        let expected_bases = vec![b'A', b'C', b'G', b'T', b'N'];
        let actual_bases: Vec<_> = [0, 1, 2, 3, 4]
            .iter()
            .map(|&b| super::u8_to_base(b, 70))
            .collect();

        assert_eq!(actual_bases, expected_bases);
    }

    #[test]
    fn extract_tiles() {
        let cbcl_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX")
            .join("Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl");
        let cbcl_header = CBCLHeader::from_path(&cbcl_path, 2).unwrap();

        let expected_bytes = vec![
            212, 254, 220, 221, 166, 108, 217, 232, 236, 221, 157, 216, 220, 220, 205, 222, 140,
            212, 157, 254, 199, 221, 237, 185, 252, 199, 237, 253, 253, 68, 237, 205, 199, 199,
            237, 109, 205, 79, 200, 220, 76, 253, 204, 253, 95, 223, 238, 78, 79, 206, 220, 152,
            220, 157, 255, 196, 207, 207, 133, 78, 236, 222, 205, 254, 237, 204, 198, 218, 236,
            204, 206, 204, 214, 207, 222, 204, 201, 221, 103, 207, 204, 196, 204, 88, 216, 205,
            222, 251, 253, 206, 206, 237, 223, 220, 205, 76, 220, 205, 232, 220,
        ];

        let mut read_buffer = vec![0u8; cbcl_header.uncompressed_size[0]];

        super::extract_tiles(&cbcl_header, 0, &mut read_buffer).unwrap();

        assert_eq!(read_buffer, expected_bytes)
    }

    #[test]
    fn extract_read_block() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();

        let expected_bq_pairs = vec![
            84, 70, 65, 70, 67, 70, 67, 70, 67, 70, 65, 70, 67, 58, 67, 70,
        ];

        let header = &novaseq_run.read_headers.get(&[1, 1]).unwrap()[0];
        let filter = &novaseq_run.filters.get(&[1, 1]).unwrap()[0];
        let pf_filter = &novaseq_run.pf_filters.get(&[1, 1]).unwrap()[0];

        let n_pf = filter.iter().map(|&b| if b { 1 } else { 0 }).sum();
        let mut bq_array = Array3::zeros((header.len(), n_pf, 2).f());
        let mut read_buffer = vec![0u8; novaseq_run.max_vec_size];
        let mut decompression_buffer = vec![0u32; 4 * novaseq_run.max_vec_size];

        super::extract_read_block(
            header,
            filter,
            pf_filter,
            &mut bq_array.slice_mut(ndarray::s![..header.len(), ..n_pf, ..]),
            &mut read_buffer,
            &mut decompression_buffer,
            0,
        );

        let bq_pairs: Vec<_> = bq_array.iter().cloned().take(16).collect();

        assert_eq!(bq_pairs, expected_bq_pairs)
    }
}
