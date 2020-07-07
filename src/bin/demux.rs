//! bcl2fastr is a program for efficient multi-threaded demultiplexing of large
//! sequencing runs (specifically from the NovaSeq instrument).

use clap::{value_t, App, Arg};
use std::path::PathBuf;
use std::str::FromStr;

use common::novaseq_run::NovaSeqRun;
use common::sample_data::read_samplesheet;
use common::write_fastq::demux_fastqs;

use log::info;
use rayon::ThreadPoolBuilder;

/// Parses command line arguments and runs demux
fn main() {
    let matches = App::new("bcl2fastr")
        .version(clap::crate_version!())
        .arg(
            Arg::with_name("run-path")
                .long("run-path")
                .help("specify path to the sequencing run folder")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("samplesheet")
                .long("samplesheet")
                .help("path to samplesheet.csv")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .long("output")
                .help("output path for fastq files")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .help("number of threads used for demultiplexing")
                .default_value("4")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("read-chunks")
                .long("read-chunks")
                .help("number of tiles to process at once while reading")
                .default_value("39")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mismatch")
                .long("mismatch")
                .help("maximum hamming distance to allow for indexes")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("compression")
                .long("compression")
                .help("compression level for gzipped output")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("verbosity")
                .short("v")
                .multiple(true)
                .help("Increase message verbosity"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .help("Silence all output"),
        )
        .arg(
            Arg::with_name("timestamp")
                .short("t")
                .help("prepend log lines with a timestamp")
                .takes_value(true)
                .possible_values(&["none", "sec", "ms", "ns"]),
        )
        .get_matches();

    let verbose = matches.occurrences_of("verbosity") as usize;
    let quiet = matches.is_present("quiet");
    let ts = matches
        .value_of("timestamp")
        .map(|v| {
            stderrlog::Timestamp::from_str(v).unwrap_or_else(|_| {
                clap::Error {
                    message: "invalid value for 'timestamp'".into(),
                    kind: clap::ErrorKind::InvalidValue,
                    info: None,
                }
                .exit()
            })
        })
        .unwrap_or(stderrlog::Timestamp::Off);

    stderrlog::new()
        .module(module_path!())
        .module("common")
        .quiet(quiet)
        .verbosity(verbose)
        .timestamp(ts)
        .init()
        .unwrap();

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

    let n_threads = value_t!(matches, "threads", usize).unwrap_or_else(|e| e.exit());
    let n_tiles = value_t!(matches, "read-chunks", usize).unwrap_or_else(|e| e.exit());
    let mismatch = value_t!(matches, "mismatch", usize).unwrap_or_else(|e| e.exit());
    let compression = value_t!(matches, "compression", u32).unwrap_or_else(|e| e.exit());

    ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build_global()
        .unwrap_or_else(|e| panic!("Error configuring global threadpool: {}", e));

    let sample_data = match read_samplesheet(samplesheet, mismatch) {
        Ok(sd) => sd,
        Err(e) => panic!("Error reading samplesheet: {}", e),
    };

    let novaseq_run = match NovaSeqRun::read_path(run_path, false) {
        Ok(n_run) => n_run,
        Err(e) => panic!("Error reading NovaSeq run: {}", e),
    };

    demux_fastqs(
        &novaseq_run,
        &sample_data,
        &output_path,
        n_tiles,
        compression,
    )
    .unwrap();

    info!("Demux complete");
}
