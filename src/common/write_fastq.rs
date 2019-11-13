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
use crate::extract_reads::{extract_read_block, extract_indices};
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
    let read_headers = novaseq_run.read_headers.get(&[lane, surface]).unwrap();
    let idx_headers = novaseq_run.index_headers.get(&[lane, surface]).unwrap();
    let filters = novaseq_run.filters.get(&[lane, surface]).unwrap();
    let pf_filters = novaseq_run.pf_filters.get(&[lane, surface]).unwrap();
    let tile_ids = novaseq_run.tile_ids.get(&[lane, surface]).unwrap();

    let mut gz_writers: Vec<_> = (1..=read_headers.len())
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

    for (i, ((filter, pf_filter), tile_id)) in filters.iter()
        .zip(pf_filters)
        .zip(tile_ids)
        .enumerate() {

        // first read the indices alone, and use them to filter the reads
        let indices = extract_indices(i, idx_headers, filter, pf_filter);

        let sample_filter: Vec<_> = indices.iter()
            .map(
                |ix|
                if let Some(_) = sample_chunk.get_sample(&ix) {
                    true
                } else {
                    false
                }
            )
            .collect();

        let samples: Vec<_> = indices.iter()
            .filter_map( |ix| sample_chunk.get_sample(&ix) )
            .collect();

        let filter: Vec<_> = filter.iter()
            .scan(0, |i, &b|
                if b {
                    *i += 1;
                    Some(sample_filter[*i - 1])
                } else {
                    Some(false)
                }
            )
            .collect();

        let pf_filter: Vec<_> = pf_filter.iter()
            .scan(0, |i, &b|
                if b {
                    *i += 1;
                    Some(sample_filter[*i - 1])
                } else {
                    Some(false)
                }
            )
            .collect();

        let locs: Vec<_> = novaseq_run.locs.iter().zip(filter.iter()).filter_map(
            |(loc, &b)| if b { Some(loc) } else { None }
        ).collect();

        let tid: Vec<_> = tile_id.iter()
            .flat_map(|(tile, n_pf)| std::iter::repeat(tile).take(*n_pf) )
            .zip(sample_filter)
            .filter_map(|(&tile, b)| if b { Some(tile) } else { None } )
            .collect();

        let n_pf = tid.len();

        for (k, (read_h, gz_writer_map)) in read_headers.iter()
            .zip(gz_writers.iter_mut())
            .enumerate() {

            let rq_array = extract_read_block(read_h, &filter, &pf_filter, i, n_pf);

            for (((rq_row, tile), loc), (sample_name, index_str)) in rq_array.axis_iter(Axis(1))
                .zip(tid.iter().cloned())
                .zip(locs.iter().cloned())
                .zip(samples.iter().cloned()) {

                let gz_writer = gz_writer_map.get_mut(&sample_name).unwrap();

                writeln!(
                    gz_writer,
                    "{}:{}:{}:{}:{} {}:N:0:{}",
                    novaseq_run.run_id,
                    lane,
                    tile,
                    loc[0],
                    loc[1],
                    k + 1,
                    index_str,
                ).unwrap();

                gz_writer.write_all(&rq_row.slice(ndarray::s![.., 0]).to_vec()).unwrap();
                gz_writer.write_all(b"\n+\n").unwrap();
                gz_writer.write_all(&rq_row.slice(ndarray::s![.., 1]).to_vec()).unwrap();
                gz_writer.write_all(b"\n").unwrap();
            }
        }
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
            if !novaseq_run.read_headers.contains_key(&[lane, surface]) {
                continue
            }

            println!("extracting lane {} surface {}", lane, surface);

            samples.par_iter()
                .panic_fuse()
                .for_each( |sample_chunk|
                    extract_sample_chunk(
                        novaseq_run, sample_chunk, output_path, lane, surface
                    )
                );
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
