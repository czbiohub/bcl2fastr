//! bcl2fastr is a program for efficient multi-threaded demultiplexing of large
//! sequencing runs (specifically from the NovaSeq instrument).

use clap::{value_t, App, Arg};
use std::path::PathBuf;
use std::str::FromStr;

use common::index_count::index_count;
use common::novaseq_run::NovaSeqRun;

use log::info;
use rayon::ThreadPoolBuilder;
use stderrlog;

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
            Arg::with_name("output")
                .long("output")
                .help("output path for index count file")
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
            Arg::with_name("top-n")
                .long("top-n")
                .help("return the top N index counts")
                .default_value("384")
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

    let output_path = PathBuf::from(matches.value_of("output").unwrap());
    if let Some(parent) = output_path.parent() {
        if !parent.is_dir() {
            panic!("Could not find output path {}", parent.display());
        }
    } else {
        panic!("output must be a file path");
    }

    let threads = value_t!(matches, "threads", usize).unwrap_or_else(|e| e.exit());
    let n_tiles = value_t!(matches, "read-chunks", usize).unwrap_or_else(|e| e.exit());
    let top_n = value_t!(matches, "top-n", usize).unwrap_or_else(|e| e.exit());

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap_or_else(|e| panic!("Error configuring global threadpool: {}", e));

    let novaseq_run = match NovaSeqRun::read_path(run_path, true) {
        Ok(n_run) => n_run,
        Err(e) => panic!("Error reading NovaSeq run: {}", e),
    };

    info!("Finished loading novaseq run");
    info!("Counting indexes");
    index_count(&novaseq_run, output_path, top_n, n_tiles).unwrap();
    info!("Done");
}
