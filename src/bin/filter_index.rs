//! index_filter is a quick script to take the output of bcl2index
//! and remove any index that matches a sample sheet

use clap::{value_t, App, Arg};
use std::{fs::File, io::prelude::*, path::PathBuf};

use common::sample_data::read_samplesheet;

/// Helper function to translate text indices to slice of ArrayView1
fn text_to_vecs(index_string: &str) -> Vec<Vec<u8>> {
    if index_string.contains("+") {
        index_string
            .split("+")
            .map(|s| s.as_bytes().to_vec())
            .collect()
    } else {
        vec![index_string.as_bytes().to_vec()]
    }
}

/// Parses command line arguments and runs demux
fn main() {
    let matches = App::new("index_filter")
        .version(clap::crate_version!())
        .arg(
            Arg::with_name("samplesheet")
                .long("samplesheet")
                .help("path to samplesheet.csv")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("input-count")
                .long("input-count")
                .help("path to unfiltered index count txt file")
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
            Arg::with_name("mismatch")
                .long("mismatch")
                .help("maximum hamming distance to allow for indexes")
                .default_value("1")
                .takes_value(true),
        )
        .get_matches();

    let samplesheet = PathBuf::from(matches.value_of("samplesheet").unwrap());
    if !samplesheet.exists() {
        panic!("Could not find samplesheet {}", samplesheet.display());
    }

    let input_path = PathBuf::from(matches.value_of("input-count").unwrap());
    if !input_path.exists() {
        panic!("Could not find input file {}", input_path.display());
    }

    let output_path = PathBuf::from(matches.value_of("output").unwrap());
    if let Some(parent) = output_path.parent() {
        if !parent.is_dir() {
            panic!("Could not find output path {}", parent.display());
        }
    } else {
        panic!("output must be a file path");
    }

    let mismatch = value_t!(matches, "mismatch", usize).unwrap_or_else(|e| e.exit());

    let sample_data = match read_samplesheet(samplesheet, mismatch) {
        Ok(sd) => sd,
        Err(e) => panic!("Error reading samplesheet: {}", e),
    };

    let mut rdr = match csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(input_path)
    {
        Ok(rdr) => rdr,
        Err(e) => panic!("Error opening index file: {}", e),
    };

    let mut out_file = match File::create(output_path) {
        Ok(out_file) => out_file,
        Err(e) => panic!("Error creating file: {}", e),
    };

    for row in rdr.records().filter_map(|r| match r {
        Ok(r) => Some(r),
        Err(e) => panic!("{}", e),
    }) {
        let indices = text_to_vecs(row.get(0).unwrap());

        if !sample_data
            .values()
            .any(|lane_samples| lane_samples.is_any_sample(&indices))
        {
            writeln!(
                &mut out_file,
                "{}\t{}",
                row.get(0).unwrap(),
                row.get(1).unwrap()
            )
            .unwrap();
        }
    }
}
