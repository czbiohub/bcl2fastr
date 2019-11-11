//! Extract only the indexes from a run and count them

use std::{
    fs::File,
    io::prelude::*,
    path::PathBuf,
};

use ndarray::Axis;
use counter::Counter;
use rayon::prelude::*;

use crate::cbcl_header_decoder::CBCLHeader;
use crate::novaseq_run::NovaSeqRun;
use crate::extract_reads::extract_reads;


fn count_tile_chunk(
    i: usize,
    headers: &[CBCLHeader],
    filter: &[bool],
    pf_filter: &[bool],
    index_ix: &[(usize, usize)],
    n_counts: usize
) -> Counter<String> {
    let reads_and_qscores = extract_reads(&headers, filter, pf_filter, i); 

    let this_count: Counter<String> = reads_and_qscores.axis_iter(Axis(1))
        .map( |rq_row| {
            index_ix.iter()
                .cloned()
                .map( |(is, ie)| rq_row.slice(ndarray::s![is..ie, 0]).to_vec() )
                .map( |v| 
                    unsafe { String::from_utf8_unchecked(v) }
                ).collect::<Vec<String>>()
                .join("+")
        }
    ).collect();

    let this_count: Counter<String> = this_count.most_common()
        .iter()
        .take(n_counts)
        .cloned()
        .collect();

    this_count
}



/// Iterate through all lanes and surfaces and count indexes
pub fn index_count(
    novaseq_run: &NovaSeqRun, output_path: PathBuf, top_n_counts: usize
) -> Result<(), &'static str> {
    let mut out_file = match File::create(output_path) {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    let top_8n_counts = top_n_counts * 8;
    let mut counts: Counter<String> = Counter::new();

    for lane in 1..=novaseq_run.run_info.flowcell_layout.lane_count {
        for surface in 1..=novaseq_run.run_info.flowcell_layout.surface_count {
            println!("indexing lane {} surface {}", lane, surface);
            let headers = novaseq_run.headers.get(&(lane, surface)).unwrap();
            let filters = novaseq_run.filters.get(&(lane, surface)).unwrap();
            let pf_filters = novaseq_run.pf_filters.get(&(lane, surface)).unwrap();

            let this_count: Counter<String> = filters.par_iter()
                .zip(pf_filters)
                .enumerate()
                .map( |(i, (filter, pf_filter))| 
                    count_tile_chunk(i, headers, filter, pf_filter, &novaseq_run.index_ix, top_8n_counts)
                ).reduce(
                    Counter::new,
                    |a, b| a + b
                );

            println!("done, adding to counts");
            counts += this_count;
        }
    }

    for (elem, freq) in counts.most_common_ordered().iter().take(top_n_counts) {
        writeln!(&mut out_file, "{}\t{}", elem, freq).unwrap();
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
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, true).unwrap();

        super::index_count(&novaseq_run, output_path, 384).unwrap()
    }

    #[test]
    #[should_panic(
        expected = r#"No such file or directory"#
    )]
    fn index_count_bad_path() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let output_path = PathBuf::from("test_data/wrong_test_output/index_counts.txt");
        let novaseq_run = NovaSeqRun::read_path(run_path, 2, true).unwrap();

        super::index_count(&novaseq_run, output_path, 384).unwrap()
    }
}
