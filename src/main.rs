//! bcl2fastr is a program for efficient multi-threaded demultiplexing of large
//! sequencing runs (specifically from the NovaSeq instrument).

use std::path::Path;
use clap::{Arg, App, value_t};

mod run_info_parser;
mod filter_decoder;
mod locs_decoder;
mod cbcl_header_decoder;
mod extract_reads;


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
        .arg(Arg::with_name("split-lanes")
             .long("split-lanes")
             .help("flag to split output by lane")).
        get_matches();

    let run_path = matches.value_of("run-path").unwrap();
    let _samplesheet = matches.value_of("samplesheet").unwrap();
    let _output = matches.value_of("output").unwrap();
    let _threads = value_t!(matches, "threads", u32).unwrap_or_else(|e| e.exit());

    if matches.is_present("split-lanes") {
        // place holder for future directives
        println!("split-lane flag is on");
    } else {
        println!("split-lane flag is off");
    }
    
    let locs_path = Path::new(run_path).join("Data/Intensities/s.locs");
    let run_info_path = Path::new(run_path).join("RunInfo.xml");
    let lane_path = Path::new(run_path).join("Data/Intensities/BaseCalls/L001");
    let tile_idces = (0, 1);    

    extract_reads::extract_reads(&locs_path, &run_info_path, &lane_path, tile_idces).unwrap();
}
