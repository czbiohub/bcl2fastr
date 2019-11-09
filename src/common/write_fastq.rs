//! Extract the reads from a run and write them out to fastq.gz files

use std::{
    fs::File,
    io::prelude::*,
    path::PathBuf,
    collections::HashMap
};

use flate2::write::GzEncoder;
use ndarray::Axis;
use rayon::prelude::*;

use crate::novaseq_run::NovaSeqRun;
use crate::extract_reads::extract_reads;
use crate::sample_data::Samples;


/// produce the correct filename format, depending on whether we are splitting lanes
fn make_filename(
    output_path: &PathBuf, sample_name: &String, lane: usize, read_num: usize
) -> PathBuf {
    if lane == 0 {
        output_path.join(format!("{}_R{}.fastq.gz", sample_name, read_num))
    } else {
        output_path.join(format!("{}_L{:03}_R{}.fastq.gz", sample_name, lane, read_num))
    }
}


/// Process a chunk of samples
fn extract_sample_chunk(
    novaseq_run: &NovaSeqRun,
    sample_chunk: &Samples,
    output_path: &PathBuf,
    lane: usize,
    surface: usize,
) {
    let headers = novaseq_run.headers.get(&(lane, surface)).unwrap();
    let filters = novaseq_run.filters.get(&(lane, surface)).unwrap();
    let pf_filters = novaseq_run.pf_filters.get(&(lane, surface)).unwrap();
    let read_ids = novaseq_run.read_ids.get(&(lane, surface)).unwrap();

    let mut gz_writers: Vec<_> = (1..=novaseq_run.read_ix.len())
        .map( |i|
            sample_chunk.sample_names.iter().map( |sample_name| {
                let out_file = match File::create(
                    make_filename(output_path, sample_name, lane, i)
                ) {
                    Ok(out_file) => out_file,
                    Err(e) => panic!("Error creating file: {}", e),
                };

                (sample_name, GzEncoder::new(out_file, flate2::Compression::new(9)))
            }).collect::<HashMap<_, _>>()
        ).collect();

    for (i, ((filter, pf_filter), rid)) in filters.iter()
        .zip(pf_filters)
        .zip(read_ids)
        .enumerate() {

        if let Ok((reads, qscores)) = extract_reads(&headers, filter, pf_filter, i) {
            for (j, sample_name, index_str, r_row, q_row) in reads.axis_iter(Axis(1))
                .zip(qscores.axis_iter(Axis(1)))
                .enumerate()
                .filter_map( |(j, (r_row, q_row))| {
                    let indices: Vec<Vec<u8>> = novaseq_run.index_ix
                        .iter()
                        .cloned()
                        .map( |(is, ie)| r_row.slice(ndarray::s![is..ie]).to_vec())
                        .collect();

                    if let Some((sample_name, index_str)) = sample_chunk.get_sample(
                        &indices
                    ) {
                        Some((j, sample_name, index_str, r_row, q_row))
                    } else {
                        None
                    }
                }) {

                for (k, (rs, re)) in novaseq_run.read_ix.iter()
                                                        .cloned()
                                                        .enumerate() {

                    let gz_writer = gz_writers[k].get_mut(&sample_name).unwrap();

                    writeln!(
                        gz_writer, "{} {}:N:0:{}", rid[j], k + 1, index_str,
                    ).unwrap();

                    gz_writer.write_all(&r_row.slice(ndarray::s![rs..re]).to_vec()).unwrap();
                    gz_writer.write_all(b"\n+\n").unwrap();
                    gz_writer.write_all(&q_row.slice(ndarray::s![rs..re]).to_vec()).unwrap();
                    gz_writer.write_all(b"\n").unwrap();
                }
            }
        };
    }
}



/// Iterate through all lanes and surfaces of a run and extract tiles in chunks
pub fn write_fastqs(
    novaseq_run: &NovaSeqRun, lane_n: usize, samples: &[Samples], output_path: &PathBuf
) -> Result<(), &'static str> {
    let lane_iter = if lane_n == 0 {
        1..=novaseq_run.run_info.flowcell_layout.lane_count
    } else {
        lane_n..=lane_n
    };

    for lane in lane_iter {
        for surface in 1..=novaseq_run.run_info.flowcell_layout.surface_count {
            if !novaseq_run.headers.contains_key(&(lane, surface)) {
                continue
            }

            println!("extracting lane {} surface {}", lane, surface);

            samples.par_iter()
                .panic_fuse()
                .map( |sample_chunk| 
                    extract_sample_chunk(
                        novaseq_run, sample_chunk, output_path, lane, surface
                    )
                )
                .for_each(drop);
        }
    }

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use crate::sample_data;

    #[test]
    fn make_filename() {
        let output_path = PathBuf::from("test_data/test_output");
        let sample_name = "sample_1".to_string();

        let file_name1 = super::make_filename(&output_path, &sample_name, 0, 1);
        assert_eq!(file_name1, output_path.join("sample_1_R1.fastq.gz"));

        let file_name2 = super::make_filename(&output_path, &sample_name, 1, 2);
        assert_eq!(file_name2, output_path.join("sample_1_L001_R2.fastq.gz"));
    }

    #[test]
    fn extract_sample_chunk() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let samplesheet_path = PathBuf::from(
            "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv"
        );
        let output_path = PathBuf::from("test_data/test_output");

        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();
        let sampledata = sample_data::read_samplesheet(samplesheet_path, 1, 1).unwrap();
        let samples = sampledata.get(&1).unwrap();

        super::extract_sample_chunk(&novaseq_run, &samples[0], &output_path, 1, 1);
    }

    #[test]
    fn write_fastqs() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let samplesheet_path = PathBuf::from(
            "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv"
        );
        let output_path = PathBuf::from("test_data/test_output");

        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();
        let sampledata = sample_data::read_samplesheet(samplesheet_path, 1, 1).unwrap();
        let samples = sampledata.get(&1).unwrap();

        super::write_fastqs(&novaseq_run, 1, samples, &output_path).unwrap();
    }
}
