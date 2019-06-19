use std::fs;
use serde_xml_rs::from_reader;


#[derive(Debug, Deserialize)]
struct Read {
    //run_path: String, // path to specific run dir
    pub number: String, // parsed from run info, outputed in fastq
    pub NumCycles: String,
    pub isIndexedRead: String,
}

#[derive(Debug, Deserialize)]
struct Reads {
    #[serde(rename = "Read", default)]
    pub Read: Vec<Read>
}

#[derive(Debug, Deserialize)]
struct Run {
    //run_path: String, // path to specific run dir
    pub Id: String, // parsed from run info, outputed in fastq
    pub Number: String,

    #[serde(rename = "Read", default)]
    pub Reads: Vec<Reads>
    
}

#[derive(Debug, Deserialize)]
struct RunInfo {
    pub Version: String,

    #[serde(rename = "Run", default)]
    pub Runs: Vec<Run>
    
}

pub fn main(run_path: String) {
    println!("reading file {}", run_path);
    let xml = fs::read_to_string(run_path).expect("error reading the file");
    let runinfo: RunInfo = from_reader(xml.as_bytes()).unwrap();
    println!("{:#?}", runinfo);
}
