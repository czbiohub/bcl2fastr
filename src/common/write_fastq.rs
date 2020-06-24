//! Extract the reads from a run and write them out to fastq.gz files

use std::{
    collections::HashMap,
    fs::{create_dir, File, OpenOptions},
    io::prelude::*,
    path::PathBuf,
};

use flate2::write::GzEncoder;
use ndarray::{Array3, ArrayView3, Axis, ShapeBuilder};
use rayon::prelude::*;

use crate::extract_reads::extract_cbcl;
use crate::novaseq_run::NovaSeqRun;
use crate::sample_data::Samples;

/// produce the correct filename format, depending on whether we are splitting lanes
fn make_filename(
    output_path: &PathBuf,
    sample_name: &String,
    sample_project: &Option<String>,
    lane: usize,
    read_num: usize,
) -> std::io::Result<PathBuf> {
    let sample_path = match sample_project {
        Some(project_name) => output_path.join(project_name),
        None => output_path.to_path_buf(),
    };

    if !sample_path.exists() {
        create_dir(&sample_path).unwrap();
    }

    if lane == 0 {
        Ok(sample_path.join(format!("{}_R{}.fastq.gz", sample_name, read_num)))
    } else {
        Ok(sample_path.join(format!(
            "{}_L{:03}_R{}.fastq.gz",
            sample_name, lane, read_num
        )))
    }
}

/// helper function to construct the report filename, depending on lane splitting
fn make_report_filename(output_path: &PathBuf, lane: usize) -> PathBuf {
    if lane == 0 {
        output_path.join(format!("barcode_report.txt"))
    } else {
        output_path.join(format!("barcode_L{:03}_report.txt", lane))
    }
}

/// remove any existing output files
fn get_sample_filepaths(
    novaseq_run: &NovaSeqRun,
    samples: &Samples,
    lane_n: usize,
    output_path: &PathBuf,
) -> std::io::Result<Vec<Vec<PathBuf>>> {
    let num_reads = novaseq_run
        .run_info
        .reads
        .iter()
        .map(|r| if r.is_indexed_read { 0 } else { 1 })
        .sum();

    let mut sample_filepaths = Vec::new();
    let mut removed_files = 0;

    for read_num in 1..=num_reads {
        let mut read_filepaths = Vec::new();

        for (sample_name, sample_project) in samples
            .sample_names
            .iter()
            .zip(samples.project_names.iter())
        {
            let file_path =
                make_filename(output_path, sample_name, sample_project, lane_n, read_num)?;

            if file_path.exists() {
                std::fs::remove_file(&file_path)?;
                removed_files += 1;
            }
            read_filepaths.push(file_path);
        }
        sample_filepaths.push(read_filepaths);
    }

    println!("removed {} files", removed_files);

    Ok(sample_filepaths)
}

/// iterate over the index arrays and count up all the reads that match sample data,
/// keeping track of whether they are an exact match to the index or not
fn count_reads(
    samples: &Samples,
    sample_i: usize,
    n_pf: usize,
    index_array: &ArrayView3<u8>,
    index_slices: &[[usize; 2]],
) -> [u64; 2] {
    let read_count = index_array
        .axis_iter(Axis(1))
        .take(n_pf)
        .par_bridge()
        .filter_map(|ix_row| {
            let indices: Vec<_> = index_slices
                .iter()
                .cloned()
                .map(|[i0, i1]| ix_row.slice(ndarray::s![i0..i1, 0]))
                .collect();

            if samples.get_sample(sample_i, &indices) {
                if !samples.is_exact(sample_i, &indices) {
                    return Some([0, 1]);
                } else {
                    return Some([1, 0]);
                }
            }

            return None;
        })
        .reduce(|| [0, 0], |a, b| [a[0] + b[0], a[1] + b[1]]);

    read_count
}

/// write the reads for a given sample to a fastq.gz file
fn write_reads(
    novaseq_run: &NovaSeqRun,
    samples: &Samples,
    sample_i: usize,
    sample_filepath: &PathBuf,
    buffer_array: &ArrayView3<u8>,
    index_array: &ArrayView3<u8>,
    index_slices: &[[usize; 2]],
    locs_vec: &[[u32; 2]],
    tile: u32,
    lane: usize,
    read_num: usize,
    compression: u32,
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

    let mut gz_writer = GzEncoder::new(out_file, flate2::Compression::new(compression));

    buffer_array
        .axis_iter(Axis(1))
        .zip(index_array.axis_iter(Axis(1)))
        .zip(locs_vec)
        .for_each(|((bq_row, ix_row), loc)| {
            let indices: Vec<_> = index_slices
                .iter()
                .cloned()
                .map(|[i0, i1]| ix_row.slice(ndarray::s![i0..i1, 0]))
                .collect();

            if samples.get_sample(sample_i, &indices) {
                write!(
                    gz_writer,
                    "{}:{}:{}:{}:{} {}:N:0:",
                    novaseq_run.run_id, lane, tile, loc[0], loc[1], read_num,
                )
                .unwrap();
                gz_writer
                    .write_all(ix_row.slice(ndarray::s![.., 0]).as_slice().unwrap())
                    .unwrap();
                gz_writer
                    .write_all(bq_row.slice(ndarray::s![.., 0]).as_slice().unwrap())
                    .unwrap();
                gz_writer.write_all(b"\n+\n").unwrap();
                gz_writer
                    .write_all(bq_row.slice(ndarray::s![.., 1]).as_slice().unwrap())
                    .unwrap();
                gz_writer.write_all(b"\n").unwrap();
            }
        });
}

/// write the read count (# total, exact, and mismatch reads) to a text file
fn write_report(
    report_filepath: &PathBuf,
    samples: &Samples,
    sample_counts: &HashMap<usize, [u64; 2]>,
) {
    let mut report_out_file = match File::create(report_filepath) {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    report_out_file
        .write_all(b"sample_name\ttotal_reads\texact_index\tindex_with_error\n")
        .unwrap();

    // sort the counts so that the output is stable
    let mut count_vec: Vec<_> = sample_counts.iter().collect();
    count_vec.sort_by_key(|(&k, _)| k);

    for (sample_i, [n_reads, m_reads]) in count_vec {
        writeln!(
            report_out_file,
            "{}\t{}\t{}\t{}",
            samples.sample_names[sample_i.clone()],
            n_reads + m_reads,
            n_reads,
            m_reads
        )
        .unwrap();
    }
}

/// Iterate through all lanes and surfaces of a run and extract tiles in chunks
pub fn demux_fastqs(
    novaseq_run: &NovaSeqRun,
    lane_n: usize,
    samples: &Samples,
    output_path: &PathBuf,
    n_chunks: usize,
    compression: u32,
) -> Result<(), &'static str> {
    // 0. check for existing files and get shared file -> path map
    let sample_files = match get_sample_filepaths(novaseq_run, samples, lane_n, output_path) {
        Ok(sample_fs) => sample_fs,
        Err(e) => panic!("Couldn't clear existing files: {}", e),
    };
    // keep track of per-sample counts and output to a report text file
    let mut sample_counts: HashMap<usize, [u64; 2]> = (0..samples.sample_names.len())
        .map(|sample_i| (sample_i, [0u64; 2]))
        .collect();
    let report_filepath = make_report_filename(output_path, lane_n);

    println!(
        "sample files: {}",
        sample_files.iter().map(|sf| sf.len()).sum::<usize>()
    );

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
    let n_cycles = novaseq_run
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
        .unwrap();

    let idx_reads: Vec<_> = novaseq_run
        .run_info
        .reads
        .iter()
        .filter(|r| r.is_indexed_read)
        .collect();

    let idx_slices: Vec<_> = idx_reads
        .iter()
        .scan(0, |k, r| {
            let t = [*k, *k + r.num_cycles];
            *k += r.num_cycles + 1;
            Some(t)
        })
        .collect();

    // leave space between index cycles for '+' and at the end for '\n'
    let n_idx_cycles = idx_reads.len() + idx_reads.iter().map(|r| r.num_cycles).sum::<usize>();

    println!("max_cycles: {}", n_cycles);

    // find the highest numbers of reads among the chunks of tiles
    let max_n_pf = novaseq_run.n_pfs.values().flatten().cloned().max().unwrap();

    println!("max_n_pf: {}", max_n_pf);

    // this array is big enough to hold all of the reads/qscores for a chunk of tiles,
    // which is basically all of the data we ever load. So as long as it fits in memory,
    // everything should be okay...
    let mut buffer_array = Array3::zeros((n_cycles, n_chunks * max_n_pf, 2).f());
    let mut index_array = Array3::zeros((n_idx_cycles, n_chunks * max_n_pf, 2).f());

    index_array
        .index_axis_mut(Axis(0), n_idx_cycles - 1)
        .fill(b'\n');

    if idx_reads.len() == 2 {
        index_array
            .index_axis_mut(Axis(0), idx_reads[0].num_cycles)
            .fill(b'+');
    }

    // preallocate vectors for loc tuples
    let mut locs_vecs = vec![Vec::with_capacity(max_n_pf); n_chunks];

    println!("buffer size: {:?}", buffer_array.raw_dim());

    for lane in lane_iter {
        for surface in novaseq_run.run_info.flowcell_layout.surface_range.clone() {
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
            let n_pfs = novaseq_run.n_pfs.get(&[lane, surface]).unwrap();

            // n_chunks defines how many tiles we extract at a time. We read all the tiles
            // in parallel within the chunk and across cycles, to maximize CPU and IO usage
            for (i, (((f_chunk, pff_chunk), tid_chunk), n_pf_chunk)) in filters
                .chunks(n_chunks)
                .zip(pf_filters.chunks(n_chunks))
                .zip(tile_ids.chunks(n_chunks))
                .zip(n_pfs.chunks(n_chunks))
                .enumerate()
            {
                println!("read chunk {}", i);

                // 1. par_iter over rows/index cycles
                // 2. per read:
                //    3. chunk_mut the array and par_iter reads into the blocks
                //    4. par_iter the reads into output files

                // beginning of this chunk
                let chunk_i = i * n_chunks;

                f_chunk
                    .par_iter()
                    .zip(&mut locs_vecs)
                    .for_each(|(filter, locs_vec)| {
                        locs_vec.clear();

                        for (loc_chunk, filt) in
                            novaseq_run.locs.chunks(2).zip(filter.iter().cloned())
                        {
                            match filt {
                                0b01 => locs_vec.push(loc_chunk[1]),
                                0b10 => locs_vec.push(loc_chunk[0]),
                                0b11 => locs_vec.extend(loc_chunk),
                                _ => (),
                            };
                        }
                    });

                println!("reading indices");
                // 1. chunk_mut the array and par_iter the indexes into it by cycle
                for (idx_vec, [idx_0, idx_1]) in idx_headers.iter().zip(idx_slices.iter().cloned())
                {
                    let mut idx_array = index_array.slice_mut(ndarray::s![idx_0..idx_1, .., ..]);

                    idx_array
                        .axis_chunks_iter_mut(Axis(1), max_n_pf)
                        .into_par_iter()
                        .zip(f_chunk)
                        .zip(pff_chunk)
                        .zip(n_pf_chunk)
                        .enumerate()
                        .for_each(|(k, (((mut ix_array, filter), pf_filter), &n_pf))| {
                            ix_array
                                .axis_iter_mut(Axis(0))
                                .into_par_iter()
                                .zip(idx_vec)
                                .for_each(|(mut byte_array, idx_h)| {
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
                                });
                        });
                }

                // 1a. count the reads for each sample
                index_array
                    .axis_chunks_iter(Axis(1), max_n_pf)
                    .zip(n_pf_chunk)
                    .for_each(|(ix_array, &n_pf)| {
                        sample_counts
                            .par_iter_mut()
                            .for_each(|(sample_i, [n_reads, m_reads])| {
                                let [lane_n, lane_m] = count_reads(
                                    samples,
                                    sample_i.clone(),
                                    n_pf,
                                    &ix_array,
                                    &idx_slices,
                                );
                                *n_reads += lane_n;
                                *m_reads += lane_m;
                            });
                    });

                // 2. per read:
                for (k, (read_h, read_files)) in read_headers.iter().zip(&sample_files).enumerate()
                {
                    println!("reading data for read {}", k + 1);
                    // 3. par_iter over cycles and read the data in
                    buffer_array
                        .axis_chunks_iter_mut(Axis(1), max_n_pf)
                        .into_par_iter()
                        .zip(f_chunk)
                        .zip(pff_chunk)
                        .zip(n_pf_chunk)
                        .enumerate()
                        .for_each(|(j, (((mut b_array, filter), pf_filter), &n_pf))| {
                            b_array
                                .axis_iter_mut(Axis(0))
                                .into_par_iter()
                                .zip(read_h)
                                .for_each(|(mut byte_array, header)| {
                                    extract_cbcl(
                                        header,
                                        if header.non_pf_clusters_excluded {
                                            &pf_filter
                                        } else {
                                            &filter
                                        },
                                        &mut byte_array.slice_mut(ndarray::s![..n_pf, ..]),
                                        chunk_i + j,
                                    )
                                });
                        });

                    println!("writing out read {}", k + 1);
                    // 4. par_iter the reads into files.
                    // It's possible/likely that n_samples >> n_threads, but they will block
                    // on i/o and so this should maximize CPU usage (maybe)
                    read_files
                        .par_iter()
                        .enumerate()
                        .for_each(|(sample_i, sample_filepath)| {
                            buffer_array
                                .axis_chunks_iter(Axis(1), max_n_pf)
                                .zip(index_array.axis_chunks_iter(Axis(1), max_n_pf))
                                .zip(&locs_vecs)
                                .zip(tid_chunk)
                                .zip(n_pf_chunk)
                                .for_each(|((((b_array, ix_array), locs_vec), &tid), &n_pf)| {
                                    write_reads(
                                        novaseq_run,
                                        samples,
                                        sample_i,
                                        sample_filepath,
                                        &b_array.slice(ndarray::s![..read_h.len(), ..n_pf, ..]),
                                        &ix_array.slice(ndarray::s![.., ..n_pf, ..]),
                                        &idx_slices,
                                        locs_vec,
                                        tid,
                                        lane,
                                        k + 1,
                                        compression,
                                    )
                                });
                        });
                }
            }
        }
    }

    write_report(&report_filepath, samples, &sample_counts);

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
        let project_name = Some("project_1".to_string());

        let file_name1 =
            super::make_filename(&output_path, &sample_name, &project_name, 0, 1).unwrap();
        assert_eq!(
            file_name1,
            output_path.join("project_1").join("sample_1_R1.fastq.gz")
        );

        let file_name2 = super::make_filename(&output_path, &sample_name, &None, 1, 2).unwrap();
        assert_eq!(file_name2, output_path.join("sample_1_L001_R2.fastq.gz"));
    }

    #[test]
    fn make_report_filename() {
        let output_path = PathBuf::from("test_data/test_output");

        let file_name1 = super::make_report_filename(&output_path, 0);
        assert_eq!(file_name1, output_path.join("barcode_report.txt"));

        let file_name2 = super::make_report_filename(&output_path, 1);
        assert_eq!(file_name2, output_path.join("barcode_L001_report.txt"));
    }

    #[test]
    fn demux_fastqs() {
        let run_path = PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX");
        let samplesheet_path =
            PathBuf::from("test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv");
        let output_path = PathBuf::from("test_data/test_output");

        let novaseq_run = NovaSeqRun::read_path(run_path, false).unwrap();
        let sampledata = sample_data::read_samplesheet(samplesheet_path, 1).unwrap();
        let samples = sampledata.get(&1).unwrap();

        super::demux_fastqs(&novaseq_run, 1, samples, &output_path, 2, 1).unwrap();
    }
}
