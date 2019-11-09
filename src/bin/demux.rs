//! bcl2fastr is a program for efficient multi-threaded demultiplexing of large
//! sequencing runs (specifically from the NovaSeq instrument).

use std::path::PathBuf;
use clap::{Arg, App, value_t};

use common::write_fastq::write_fastqs;
use common::novaseq_run::NovaSeqRun;
use common::sample_data::read_samplesheet;

use rayon::ThreadPoolBuilder;


/// Parses command line arguments and runs demux
fn main() {
    let matches = App::new("bcl2fastr")
        .version(clap::crate_version!())
        .arg(Arg::with_name("run-path")
            .long("run-path")
            .help("specify path to the sequencing run folder")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("samplesheet")
            .long("samplesheet")
            .help("path to samplesheet.csv")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("output")
            .long("output")
            .help("output path for fastq files")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("threads")
            .long("threads")
            .help("number of threads used for demultiplexing")
            .default_value("4")
            .takes_value(true))
        .arg(Arg::with_name("tile-chunk")
            .long("tile-chunk")
            .help("number of tiles to extract at a time")
            .default_value("2")
            .takes_value(true))
        .arg(Arg::with_name("mismatch")
            .long("mismatch")
            .help("maximum hamming distance to allow for indexes")
            .default_value("1")
            .takes_value(true))
        .get_matches();

    let run_path = PathBuf::from(matches.value_of("run-path").unwrap());
    if !run_path.exists() {
        panic!("Could not find run path {}", run_path.display());
    }

    let samplesheet = PathBuf::from(matches.value_of("samplesheet").unwrap());
    if !samplesheet.exists() {
        panic!("Could not find samplesheet {}", samplesheet.display());
    }

    let output_path = PathBuf::from(matches.value_of("output").unwrap());
    if !output_path.exists() {
        panic!("Could not find output path {}", output_path.display());
    }
    if !output_path.is_dir() {
        panic!("Output path {} is not a directory", output_path.display());
    }

    let threads = value_t!(matches, "threads", usize).unwrap_or_else(|e| e.exit());
    let tile_chunk = value_t!(matches, "tile-chunk", usize).unwrap_or_else(|e| e.exit());
    let mismatch = value_t!(matches, "mismatch", usize).unwrap_or_else(|e| e.exit());

    ThreadPoolBuilder::new().num_threads(threads).build_global()
        .unwrap_or_else(|e| panic!("Error configuring global threadpool: {}", e));

    let sample_data = match read_samplesheet(samplesheet, mismatch, threads) {
        Ok(sd) => sd,
        Err(e) => panic!("Error reading samplesheet: {}", e),
    };

    let novaseq_run = match NovaSeqRun::read_path(run_path, tile_chunk, false) {
        Ok(n_run) => n_run,
        Err(e) => panic!("Error reading NovaSeq run: {}", e),
    };

    for (lane, sample_vec) in sample_data {
        write_fastqs(&novaseq_run, lane.clone(), &sample_vec, &output_path).unwrap();
    }
}
