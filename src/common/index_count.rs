//! Extract only the indexes from a run and count them

use std::{cmp::Ordering, fs::File, io::prelude::*, path::PathBuf};

use counter::Counter;
use rayon::prelude::*;

use ndarray::{Array3, Axis, ShapeBuilder};

use crate::cbcl_header_decoder::CBCLHeader;
use crate::extract_reads::extract_cbcl;
use crate::novaseq_run::NovaSeqRun;

use log::{debug, info};

fn count_tile_chunk(
    tile_i: usize,
    headers: &[Vec<CBCLHeader>],
    filter: &[u8],
    pf_filter: &[u8],
    n_pf: usize,
    n_counts: usize,
) -> Counter<Vec<u8>> {
    let n_idx_cycles: usize = headers.iter().map(|h| h.len()).sum();

    let mut index_array = Array3::zeros((n_idx_cycles + headers.len() - 1, n_pf, 2).f());
    if headers.len() == 2 {
        index_array
            .index_axis_mut(Axis(0), headers[0].len())
            .fill(b'+');
    }

    let mut j = 0;
    for idx_vec in headers {
        let mut idx_array = index_array.slice_mut(ndarray::s![j..j + idx_vec.len(), ..n_pf, ..]);
        j += idx_vec.len() + 1;

        for (mut byte_array, idx_h) in idx_array.axis_iter_mut(Axis(0)).zip(idx_vec) {
            extract_cbcl(
                idx_h,
                if idx_h.non_pf_clusters_excluded {
                    &pf_filter
                } else {
                    &filter
                },
                &mut byte_array.slice_mut(ndarray::s![..n_pf, ..]),
                tile_i,
            )
        }
    }

    let this_count: Counter<Vec<u8>> = index_array
        .index_axis(Axis(2), 0)
        .axis_iter(Axis(1))
        .map(|ix| ix.to_vec())
        .collect();

    this_count
        .most_common()
        .iter()
        .take(n_counts)
        .cloned()
        .collect()
}

/// Iterate through all lanes and surfaces and count indexes
pub fn index_count(
    novaseq_run: &NovaSeqRun,
    output_path: PathBuf,
    top_n_counts: usize,
    k_fold: usize,
) -> Result<(), &'static str> {
    info!("writing to {}", output_path.display());
    let mut out_file = match File::create(output_path) {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    let top_kn_counts = top_n_counts * k_fold;
    let mut counters = Vec::new();

    info!("Counting indexes");
    for lane in 1..=novaseq_run.run_info.flowcell_layout.lane_count {
        debug!("Starting lane {}", lane);
        for surface in novaseq_run.run_info.flowcell_layout.surface_range.clone() {
            debug!("Starting surface {}", surface);

            let filters = novaseq_run.filters.get(&[lane, surface]).unwrap();
            let pf_filters = novaseq_run.pf_filters.get(&[lane, surface]).unwrap();
            let idx_headers = novaseq_run.index_headers.get(&[lane, surface]).unwrap();
            let n_pfs = novaseq_run.n_pfs.get(&[lane, surface]).unwrap();

            let this_count: Counter<Vec<u8>> = filters
                .par_iter()
                .zip(pf_filters)
                .zip(n_pfs)
                .enumerate()
                .map(|(i, ((filter, pf_filter), &n_pf))| {
                    count_tile_chunk(i, idx_headers, filter, pf_filter, n_pf, top_kn_counts)
                })
                .reduce(Counter::new, |a, b| a + b);

            debug!("Done with {} - {}, adding to counts", lane, surface);
            counters.push(this_count);
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

        super::index_count(&novaseq_run, output_path, 384, 4).unwrap()
    }

    #[test]
    #[should_panic(expected = r#"No such file or directory"#)]
    fn index_count_bad_path() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let output_path = PathBuf::from("test_data/wrong_test_output/index_counts.txt");
        let novaseq_run = NovaSeqRun::read_path(run_path, true).unwrap();

        super::index_count(&novaseq_run, output_path, 384, 4).unwrap()
    }
}
