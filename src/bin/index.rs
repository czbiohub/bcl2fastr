//! bcl2fastr is a program for efficient multi-threaded demultiplexing of large
//! sequencing runs (specifically from the NovaSeq instrument).

use std::path::PathBuf;
use clap::{Arg, App, value_t};

use common::index_count::index_count;
use common::novaseq_run::NovaSeqRun;

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
        .arg(Arg::with_name("top-n")
            .long("top-n")
            .help("return the top N index counts")
            .default_value("384")
            .takes_value(true))
        .get_matches();

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
    let tile_chunk = value_t!(matches, "tile-chunk", usize).unwrap_or_else(|e| e.exit());
    let top_n = value_t!(matches, "top-n", usize).unwrap_or_else(|e| e.exit());

    ThreadPoolBuilder::new().num_threads(threads).build_global()
        .unwrap_or_else(|e| panic!("Error configuring global threadpool: {}", e));

    let novaseq_run = match NovaSeqRun::read_path(run_path, tile_chunk, true) {
        Ok(n_run) => n_run,
        Err(e) => panic!("Error reading NovaSeq run: {}", e),
    };

    println!("Finished loading novaseq run");
    println!("Counting indexes");
    index_count(&novaseq_run, output_path, top_n).unwrap();
}
