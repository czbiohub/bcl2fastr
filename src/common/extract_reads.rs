//! Extract and decompress a set of tiles from a vector of cbcl files.

use std::{fs::File, io::prelude::*, io::SeekFrom, ptr::write};

use flate2::read::GzDecoder;
use ndarray::{ArrayViewMut2, Axis};

use crate::base_decoder::{B_MAP_01, B_MAP_10, Q_MAP_01, Q_MAP_10};
use crate::cbcl_header_decoder::CBCLHeader;

/// extract multiple tiles from a CBCL file and write them into the array
fn extract_tiles(
    header: &CBCLHeader,
    tile_i: usize,
    bq_cycle: &mut ArrayViewMut2<u8>,
    filter: &[u8],
) -> std::io::Result<()> {
    let start_pos = header.start_pos[tile_i];
    let compressed_size = header.compressed_size[tile_i];

    // open file and seek to start position
    let mut cbcl = File::open(&header.cbcl_path)?;
    cbcl.seek(SeekFrom::Start(start_pos))?;

    let n_cycles = bq_cycle.strides()[0] as usize;
    let mut read_slice = bq_cycle.slice_mut(ndarray::s![.., 0]).as_mut_ptr();
    let mut qscore_slice = bq_cycle.slice_mut(ndarray::s![.., 1]).as_mut_ptr();

    // use GzDecoder to decompress
    let gz = GzDecoder::new(cbcl.take(compressed_size));
    for (byte, f) in gz.bytes().zip(filter) {
        let c = match byte {
            Ok(c) => c as usize,
            Err(e) => return Err(e),
        };
        match f {
            3 => unsafe {
                write(read_slice, B_MAP_10[c]);
                read_slice = read_slice.add(n_cycles);
                write(read_slice, B_MAP_01[c]);
                read_slice = read_slice.add(n_cycles);
                write(qscore_slice, Q_MAP_10[c]);
                qscore_slice = qscore_slice.add(n_cycles);
                write(qscore_slice, Q_MAP_01[c]);
                qscore_slice = qscore_slice.add(n_cycles);
            },
            1 => unsafe {
                write(read_slice, B_MAP_01[c]);
                read_slice = read_slice.add(n_cycles);
                write(qscore_slice, Q_MAP_01[c]);
                qscore_slice = qscore_slice.add(n_cycles);
            },
            2 => unsafe {
                write(read_slice, B_MAP_10[c]);
                read_slice = read_slice.add(n_cycles);
                write(qscore_slice, Q_MAP_10[c]);
                qscore_slice = qscore_slice.add(n_cycles);
            },
            _ => (),
        }
    }

    Ok(())
}

/// just read a lot of data into one cycle
pub fn extract_cbcl(
    header: &CBCLHeader,
    filter: &[u8],
    bq_cycle: &mut ArrayViewMut2<u8>,
    tile_i: usize,
) {
    if let Ok(_) = extract_tiles(header, tile_i, bq_cycle, filter) {
        return;
    } else {
        bq_cycle.index_axis_mut(Axis(1), 0).fill(b'N');
        bq_cycle.index_axis_mut(Axis(1), 1).fill(b'#');
    }
}

#[cfg(test)]
mod tests {
    use ndarray::{Array2, Array3, Axis, ShapeBuilder};
    use std::path::PathBuf;

    use crate::cbcl_header_decoder::CBCLHeader;
    use crate::filter_decoder::filter_decoder;
    use crate::novaseq_run::NovaSeqRun;

    #[test]
    fn extract_tiles() {
        let cbcl_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX")
            .join("Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl");
        let cbcl_header = CBCLHeader::from_path(&cbcl_path).unwrap();
        let filter_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX")
            .join("Data/Intensities/BaseCalls/L001/s_1_1101.filter");
        let cbcl_filter = filter_decoder(&filter_path).unwrap();
        let n_pf = cbcl_filter
            .iter()
            .map(|&b| [0, 1, 1, 2][b as usize])
            .sum::<usize>();

        let expected_bytes = vec![
            84, 70, 65, 70, 67, 70, 67, 70, 67, 70, 65, 70, 67, 58, 67, 70, 65, 58, 71, 70, 65, 70,
            71, 70, 67, 70, 67, 70, 67, 70, 65, 58, 67, 70, 65, 70, 67, 70, 65, 70, 67, 70, 67, 70,
            65, 70, 71, 70, 67, 70, 65, 70, 65, 58, 65, 44, 67, 70, 67, 70, 67, 58, 71, 70, 84, 70,
            65, 70, 67, 70, 67, 70, 67, 70, 71, 70, 67, 58, 84, 58, 65, 70, 84, 70, 65, 70, 67, 70,
            71, 70, 67, 70, 84, 70, 67, 70, 84, 70, 65, 44, 67, 70, 71, 70, 67, 70, 65, 70, 65, 70,
            84, 44, 65, 70, 67, 70, 71, 70, 67, 70, 67, 70, 65, 70, 84, 70, 65, 70, 65, 70, 67, 70,
            65, 70, 67, 70, 84, 70, 65, 70, 65, 70, 67, 70, 84, 70, 84, 70, 84, 70, 67, 70, 71, 70,
            65, 44, 84, 70,
        ];

        let mut read_array: Array2<u8> = Array2::zeros((n_pf, 2));

        super::extract_tiles(&cbcl_header, 0, &mut read_array.view_mut(), &cbcl_filter).unwrap();

        assert_eq!(read_array.into_raw_vec(), expected_bytes);
    }

    #[test]
    fn extract_read_block() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let novaseq_run = NovaSeqRun::read_path(run_path, false).unwrap();

        let expected_bq_pairs = vec![
            vec![
                84, 70, 65, 70, 67, 70, 67, 70, 67, 70, 65, 70, 67, 58, 67, 70,
            ],
            vec![
                65, 70, 71, 70, 65, 70, 84, 70, 84, 70, 84, 70, 67, 44, 84, 70,
            ],
            vec![
                84, 70, 65, 70, 67, 70, 67, 70, 67, 70, 65, 70, 65, 70, 65, 70,
            ],
            vec![
                65, 70, 67, 70, 84, 70, 84, 70, 65, 70, 84, 70, 65, 70, 67, 70,
            ],
        ];

        let headers = &novaseq_run.read_headers.get(&[1, 1]).unwrap()[0];
        let filter = &novaseq_run.filters.get(&[1, 1]).unwrap()[0];
        let pf_filter = &novaseq_run.pf_filters.get(&[1, 1]).unwrap()[0];

        let n_pf = filter.iter().map(|&b| [0, 1, 1, 2][b as usize]).sum();
        let mut bq_array = Array3::zeros((headers.len(), n_pf, 2).f());

        for ((mut byte_array, read_h), exp_bq) in bq_array
            .axis_iter_mut(Axis(0))
            .zip(headers)
            .zip(expected_bq_pairs)
        {
            super::extract_cbcl(
                read_h,
                if read_h.non_pf_clusters_excluded {
                    &pf_filter
                } else {
                    &filter
                },
                &mut byte_array,
                0,
            );
            let bq_pairs: Vec<_> = byte_array.iter().cloned().take(16).collect();
            assert_eq!(bq_pairs, exp_bq);
        }
    }
}
