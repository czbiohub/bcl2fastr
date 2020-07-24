//! Extract only the indexes from a run and count them

use std::{cmp::Ordering, fs::File, io::prelude::*, path::PathBuf};

use counter::Counter;
use rayon::prelude::*;

use ndarray::{Array3, Axis, ShapeBuilder};

use crate::extract_reads::extract_cbcl;
use crate::novaseq_run::NovaSeqRun;

use log::{debug, info};

/// Iterate through all lanes and surfaces and count indexes
pub fn index_count(
    novaseq_run: &NovaSeqRun,
    output_path: PathBuf,
    top_n_counts: usize,
    n_tiles: usize,
) -> Result<(), &'static str> {
    info!("writing to {}", output_path.display());
    let mut out_file = match File::create(output_path) {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    let idx_reads: Vec<_> = novaseq_run
        .run_info
        .reads
        .iter()
        .filter(|r| r.is_indexed_read)
        .collect();

    // slices into the index array, including padding for extra characters needed for output
    let idx_slices: Vec<_> = idx_reads
        .iter()
        .scan(0, |k, r| {
            let t = [*k, *k + r.num_cycles];
            *k += r.num_cycles + 1;
            Some(t)
        })
        .collect();

    // leave space between index cycles for '+'
    let n_idx_cycles = idx_reads.len() + idx_reads.iter().map(|r| r.num_cycles).sum::<usize>() - 1;

    // find the highest numbers of reads among the chunks of tiles
    let max_n_pf = novaseq_run.n_pfs.values().flatten().cloned().max().unwrap();

    debug!("max_n_pf: {}", max_n_pf);

    // create an array that can hold all of the indexes for a chunk of tiles
    let mut index_array = Array3::zeros((n_idx_cycles, n_tiles * max_n_pf, 2).f());

    // fill the row between index reads with '+' for output
    if idx_reads.len() == 2 {
        index_array
            .index_axis_mut(Axis(0), idx_reads[0].num_cycles)
            .fill(b'+');
    }

    let mut counters = Vec::new();

    info!("Counting indexes");
    for lane in 1..=novaseq_run.run_info.flowcell_layout.lane_count {
        for surface in novaseq_run.run_info.flowcell_layout.surface_range.clone() {
            debug!("Counting indexes for lane {} surface {}", lane, surface);

            let idx_headers = novaseq_run.index_headers.get(&[lane, surface]).unwrap();
            let filters = novaseq_run.filters.get(&[lane, surface]).unwrap();
            let pf_filters = novaseq_run.pf_filters.get(&[lane, surface]).unwrap();
            let n_pfs = novaseq_run.n_pfs.get(&[lane, surface]).unwrap();

            // n_tiles defines how many tiles we extract at a time. We read across all cycles
            for (i, ((f_chunk, pff_chunk), n_pf_chunk)) in filters
                .chunks(n_tiles)
                .zip(pf_filters.chunks(n_tiles))
                .zip(n_pfs.chunks(n_tiles))
                .enumerate()
            {
                debug!("Read chunk {}", i);

                // 1. par_iter over rows/index cycles
                // 2. count the indexes and store the counter

                // beginning of this chunk
                let chunk_i = i * n_tiles;

                debug!("Reading indices");
                // 1. chunk_mut the array and par_iter the indexes into it by cycle
                for (idx_vec, [idx_0, idx_1]) in idx_headers.iter().zip(idx_slices.iter().cloned())
                {
                    let mut idx_array = index_array.slice_mut(ndarray::s![idx_0..idx_1, .., ..]);

                    // outer loop is by cycle, to spread out file IO
                    idx_array
                        .axis_iter_mut(Axis(0))
                        .into_par_iter()
                        .zip(idx_vec)
                        .for_each(|(mut cycle_array, idx_h)| {
                            // inner loop is by tile (axis 1 in full ndarray)
                            cycle_array
                                .axis_chunks_iter_mut(Axis(0), max_n_pf)
                                .zip(f_chunk)
                                .zip(pff_chunk)
                                .zip(n_pf_chunk)
                                .enumerate()
                                .for_each(|(k, (((mut byte_array, filter), pf_filter), &n_pf))| {
                                    extract_cbcl(
                                        idx_h,
                                        if idx_h.non_pf_clusters_excluded {
                                            &pf_filter
                                        } else {
                                            &filter
                                        },
                                        &mut byte_array.slice_mut(ndarray::s![..n_pf, ..]),
                                        chunk_i + k,
                                    );
                                })
                        });
                }

                let this_count: Counter<Vec<u8>> = index_array
                    .axis_chunks_iter(Axis(1), max_n_pf)
                    .into_par_iter()
                    .zip(n_pf_chunk)
                    .map(|(ix_chunk, &n_pf)| {
                        ix_chunk
                            .slice(ndarray::s![.., ..n_pf, 0])
                            .axis_iter(Axis(1))
                            .map(|ix_row| ix_row.to_vec())
                            .collect::<Counter<Vec<u8>>>()
                    })
                    .reduce(Counter::new, |a, b| a + b);

                debug!("Done with {} - {}, adding to counts", lane, surface);
                counters.push(this_count);
            }
        }
        debug!("Lane {} complete", lane);
    }

    debug!("Merging counts");
    let counts = counters
        .par_iter()
        .cloned()
        .reduce(Counter::new, |a, b| a + b);

    debug!("Sorting");
    let mut count_vec: Vec<_> = counts.iter().collect();
    count_vec.par_sort_by(|(e1, f1), (e2, f2)| match f2.cmp(f1) {
        Ordering::Equal => e1.cmp(e2),
        unequal => unequal,
    });

    debug!("Writing");
    for (elem, freq) in count_vec.iter().take(top_n_counts) {
        writeln!(
            &mut out_file,
            "{}\t{}",
            String::from_utf8(elem.to_vec()).unwrap(),
            freq
        )
        .unwrap();
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn index_count() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let output_path = PathBuf::from("test_data/test_output/index_counts.txt");
        let novaseq_run = NovaSeqRun::read_path(run_path, true).unwrap();

        super::index_count(&novaseq_run, output_path, 384, 2).unwrap()
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn index_count_bad_path() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let output_path = PathBuf::from("test_data/wrong_test_output/index_counts.txt");
        let novaseq_run = NovaSeqRun::read_path(run_path, true).unwrap();

        super::index_count(&novaseq_run, output_path, 384, 2).unwrap()
    }
}
