//! Extract the reads from a run and write them out to fastq.gz files

use std::{fs::OpenOptions, io::prelude::*, path::PathBuf};

use flate2::write::GzEncoder;
use ndarray::{Array3, ArrayView3, Axis, ShapeBuilder};
use rayon::prelude::*;

use crate::cbcl_header_decoder::CBCLHeader;
use crate::extract_reads::extract_read_block;
use crate::novaseq_run::NovaSeqRun;
use crate::sample_data::Samples;

/// produce the correct filename format, depending on whether we are splitting lanes
fn make_filename(
    output_path: &PathBuf,
    sample_name: &String,
    lane: usize,
    read_num: usize,
) -> PathBuf {
    if lane == 0 {
        output_path.join(format!("{}_R{}.fastq.gz", sample_name, read_num))
    } else {
        output_path.join(format!(
            "{}_L{:03}_R{}.fastq.gz",
            sample_name, lane, read_num
        ))
    }
}

/// remove any existing output files
fn get_sample_filepaths<'a>(
    novaseq_run: &NovaSeqRun,
    samples: &'a Samples,
    lane_n: usize,
    output_path: &PathBuf,
) -> std::io::Result<Vec<PathBuf>> {
    let num_reads = novaseq_run
        .run_info
        .reads
        .iter()
        .map(|r| if r.is_indexed_read { 0 } else { 1 })
        .sum();

    let mut sample_filepaths = Vec::new();
    let mut removed_files = 0;

    for read_num in 1..=num_reads {
        for sample_name in samples.sample_names.iter() {
            let file_path = make_filename(output_path, sample_name, lane_n, read_num);
            if file_path.exists() {
                std::fs::remove_file(&file_path)?;
                removed_files += 1;
            }
            sample_filepaths.push(file_path);
        }
    }

    println!("removed {} files", removed_files);

    Ok(sample_filepaths)
}

fn get_samples(
    samples: &Samples,
    bq_array: &ArrayView3<u8>,
    idx_headers: &[Vec<CBCLHeader>],
    n_pf: usize,
) -> Vec<Option<[usize; 3]>> {
    // find which indices correspond to samples in any chunk
    let sample_options: Vec<_> = (0..n_pf)
        .map(|j| {
            samples.get_sample(
                &idx_headers
                    .iter()
                    .scan(0, |k, idx_h| {
                        *k += idx_h.len();
                        Some(bq_array.slice(ndarray::s![*k - idx_h.len()..*k, j, 0]))
                    })
                    .collect::<Vec<_>>()[..],
            )
        })
        .collect();

    sample_options
}

fn write_reads<'a>(
    novaseq_run: &NovaSeqRun,
    samples: &Samples,
    sample_i: usize,
    sample_filepath: &PathBuf,
    sample_options: &[Vec<Option<[usize; 3]>>],
    buffer_array: &ArrayView3<u8>,
    locs_chunk: &[Vec<[u32; 2]>],
    tid_chunk: &[Vec<u32>],
    max_n_pf: usize,
    lane: usize,
    read_num: usize,
) {
    // create gz writer for this sample, or open for appending
    let out_file = match OpenOptions::new()
        .create(true)
        .append(true)
        .open(sample_filepath)
    {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    let mut gz_writer = GzEncoder::new(out_file, flate2::Compression::new(4));

    // could par_iter here? seems like overkill
    buffer_array
        .axis_chunks_iter(Axis(1), max_n_pf)
        .zip(locs_chunk)
        .zip(tid_chunk)
        .zip(sample_options)
        .for_each(|(((bq_array, locs), tid), sample_opt)| {
            // zip will stop when all samples are used, so don't need to slice bq_array
            for (rq_row, tile, loc, idx_i, idx_i2) in bq_array
                .axis_iter(Axis(1))
                .zip(tid.iter().cloned())
                .zip(locs.iter().cloned())
                .zip(sample_opt.iter().cloned())
                .filter_map(|(((rq_row, tile), loc), sample_o)| {
                    if let Some([this_sample, idx_i, idx_i2]) = sample_o {
                        if this_sample == sample_i {
                            Some((rq_row, tile, loc, idx_i, idx_i2))
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
            {
                writeln!(
                    gz_writer,
                    "{}:{}:{}:{}:{} {}:N:0:{}{}",
                    novaseq_run.run_id,
                    lane,
                    tile,
                    loc[0],
                    loc[1],
                    read_num,
                    samples.index[idx_i],
                    samples.index2[idx_i2]
                )
                .unwrap();

                gz_writer
                    .write_all(rq_row.slice(ndarray::s![.., 0]).as_slice().unwrap())
                    .unwrap();
                gz_writer.write_all(b"\n+\n").unwrap();
                gz_writer
                    .write_all(rq_row.slice(ndarray::s![.., 1]).as_slice().unwrap())
                    .unwrap();
                gz_writer.write_all(b"\n").unwrap();
            }
        });
}

/// Iterate through all lanes and surfaces of a run and extract tiles in chunks
pub fn demux_fastqs(
    novaseq_run: &NovaSeqRun,
    lane_n: usize,
    samples: &Samples,
    output_path: &PathBuf,
    n_chunks: usize,
) -> Result<(), &'static str> {
    // 0. check for existing files and get shared file -> path map
    let sample_files = match get_sample_filepaths(novaseq_run, samples, lane_n, output_path) {
        Ok(sample_fs) => sample_fs,
        Err(e) => panic!("Couldn't clear existing files: {}", e),
    };

    println!("sample files: {}", sample_files.len());

    // if split-lanes is true, we will get a single lane number. Otherwise, lane_n is 0
    // and we should process all the lanes
    let lane_iter = if lane_n == 0 {
        1..=novaseq_run.run_info.flowcell_layout.lane_count
    } else {
        lane_n..=lane_n
    };

    // compute max array depth needed for data. There is a chance that the total index
    // length is longer than the longest read (for instance, in test data) so we account
    // for that here
    let max_cycles = std::cmp::max(
        novaseq_run
            .run_info
            .reads
            .iter()
            .filter_map(|r| {
                if r.is_indexed_read {
                    None
                } else {
                    Some(r.num_cycles)
                }
            })
            .max()
            .unwrap(),
        novaseq_run
            .run_info
            .reads
            .iter()
            .filter_map(|r| {
                if r.is_indexed_read {
                    Some(r.num_cycles)
                } else {
                    None
                }
            })
            .sum(),
    );

    println!("max_cycles: {}", max_cycles);

    // find the highest numbers of reads among the chunks of tiles
    let max_n_pf = novaseq_run
        .tile_ids
        .values()
        .flat_map(|ls_tile_ids| {
            ls_tile_ids.chunks(n_chunks).flat_map(|tids| {
                tids.iter()
                    .map(|tid| tid.iter().map(|(_, n_pf)| n_pf).sum())
            })
        })
        .max()
        .unwrap();

    println!("max_n_pf: {}", max_n_pf);

    // this array is big enough to hold all of the reads/qscores for a chunk of tiles,
    // which is basically all of the data we ever load. So as long as it fits in memory,
    // everything should be okay...
    let mut buffer_array = Array3::zeros((max_cycles, n_chunks * max_n_pf, 2).f());

    // pre-allocating arrays and vectors to re-use during reading
    let mut read_buffers = vec![vec![0u8; novaseq_run.max_vec_size]; n_chunks];
    let mut decomp_buffers = vec![vec![0u32; 4 * novaseq_run.max_vec_size]; n_chunks];

    println!("buffer size: {:?}", buffer_array.raw_dim());

    for lane in lane_iter {
        for surface in 1..=novaseq_run.run_info.flowcell_layout.surface_count {
            // check to make sure the data is here. Only relevant for testing
            if !novaseq_run.read_headers.contains_key(&[lane, surface]) {
                continue;
            }

            println!("extracting lane {} surface {}", lane, surface);

            let read_headers = novaseq_run.read_headers.get(&[lane, surface]).unwrap();
            let idx_headers = novaseq_run.index_headers.get(&[lane, surface]).unwrap();
            let filters = novaseq_run.filters.get(&[lane, surface]).unwrap();
            let pf_filters = novaseq_run.pf_filters.get(&[lane, surface]).unwrap();
            let tile_ids = novaseq_run.tile_ids.get(&[lane, surface]).unwrap();

            // n_chunks defines how many tiles we extract at a time, where tiles are
            // already merged into meta-tiles. We sequentially read these chunks,
            // processing each one in parallel across all the threads
            for (i, ((f_chunk, pff_chunk), tid_tuples)) in filters
                .chunks(n_chunks)
                .zip(pf_filters.chunks(n_chunks))
                .zip(tile_ids.chunks(n_chunks))
                .enumerate()
            {
                println!("read chunk {}", i);

                // 1. chunk_mut the index arrays and par_iter indexes in
                // 2. par_iter over data chunks to get option vecs from indexes
                // 3. per read:
                //    4. chunk_mut the array and par_iter reads into the blocks
                //    5. par_iter the reads into output files

                // beginning of this chunk
                let chunk_i = i * n_chunks;

                let locs_chunk: Vec<_> = f_chunk
                    .par_iter()
                    .map(|filter| {
                        novaseq_run
                            .locs
                            .iter()
                            .zip(filter.iter())
                            .filter_map(|(loc, &b)| if b { Some(loc.clone()) } else { None })
                            .collect::<Vec<_>>()
                    })
                    .collect();

                let tid_chunk: Vec<_> = tid_tuples
                    .par_iter()
                    .map(|tid| {
                        tid.iter()
                            .flat_map(|(tile, n_pf)| std::iter::repeat(tile.clone()).take(*n_pf))
                            .collect::<Vec<_>>()
                    })
                    .collect();

                let n_pfs: Vec<_> = tid_chunk.iter().map(|tid| tid.len()).collect();

                println!("reading indices");
                // 1. chunk_mut the array and par_iter the indexes in
                buffer_array
                    .axis_chunks_iter_mut(Axis(1), max_n_pf)
                    .into_par_iter()
                    .zip(&mut read_buffers)
                    .zip(&mut decomp_buffers)
                    .zip(f_chunk)
                    .zip(pff_chunk)
                    .zip(&n_pfs)
                    .enumerate()
                    .for_each(
                        // b_a: byte_array, r_b: read_buffer, d_b: decompression_buffer
                        |(j, (((((mut b_a, r_b), d_b), filter), pf_filter), &n_pf))| {
                            let mut k = 0;
                            for idx_h in idx_headers {
                                extract_read_block(
                                    idx_h,
                                    &filter,
                                    &pf_filter,
                                    &mut b_a.slice_mut(ndarray::s![k..k + idx_h.len(), ..n_pf, ..]),
                                    r_b,
                                    d_b,
                                    chunk_i + j,
                                );
                                k += idx_h.len();
                            }
                        },
                    );

                println!("parsing {} samples", samples.sample_names.len());

                // make one vector of options and share it out
                let mut sample_options = Vec::with_capacity(n_chunks);

                // 2. par_iter over the data chunks to get option vecs from indexes
                buffer_array
                    .axis_chunks_iter(Axis(1), max_n_pf)
                    .into_par_iter()
                    .zip(&n_pfs)
                    .map(|(bq_array, n_pf)| {
                        get_samples(&samples, &bq_array, idx_headers, n_pf.clone())
                    })
                    .collect_into_vec(&mut sample_options);

                // 3. per read:
                for (k, read_h) in read_headers.iter().enumerate() {
                    println!("reading data for read {}", k + 1);
                    // 4. chunk_mut the array and par_iter reads in
                    buffer_array
                        .axis_chunks_iter_mut(Axis(1), max_n_pf)
                        .into_par_iter()
                        .zip(&mut read_buffers)
                        .zip(&mut decomp_buffers)
                        .zip(f_chunk)
                        .zip(pff_chunk)
                        .zip(&n_pfs)
                        .enumerate()
                        .for_each(
                            // b_a: byte_array, r_b: read_buffer, d_b: decomp_buffer
                            |(j, (((((mut b_a, r_b), d_b), filter), pf_filter), &n_pf))| {
                                extract_read_block(
                                    read_h,
                                    &filter,
                                    &pf_filter,
                                    &mut b_a.slice_mut(ndarray::s![..read_h.len(), ..n_pf, ..]),
                                    r_b,
                                    d_b,
                                    chunk_i + j,
                                )
                            },
                        );

                    println!("writing out read {}", k + 1);
                    // 5. par_iter the reads into files.
                    // It's possible that n_samples >> n_threads, but they will block
                    // on i/o and so this should maximize CPU usage (maybe)
                    sample_files
                        .par_iter()
                        .enumerate()
                        .for_each(|(sample_i, sample_filepath)| {
                            write_reads(
                                novaseq_run,
                                samples,
                                sample_i,
                                sample_filepath,
                                &sample_options,
                                &buffer_array.slice(ndarray::s![..read_h.len(), .., ..]),
                                &locs_chunk,
                                &tid_chunk,
                                max_n_pf,
                                lane,
                                k + 1,
                            )
                        });
                }
            }
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
    fn demux_fastqs() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let samplesheet_path =
            PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv");
        let output_path = PathBuf::from("test_data/test_output");

        let novaseq_run = NovaSeqRun::read_path(run_path, 2, false).unwrap();
        let sampledata = sample_data::read_samplesheet(samplesheet_path, 1).unwrap();
        let samples = sampledata.get(&1).unwrap();

        super::demux_fastqs(&novaseq_run, 1, samples, &output_path, 2).unwrap();
    }
}
