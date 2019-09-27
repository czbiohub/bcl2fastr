//! Deserializes the `RunInfo.xml` file from a NovaSeq run into a useful struct
//! of information about the sequencing run.

use std::{
    fs,
    path::Path,
};
use serde::{de, Deserialize};
use serde_xml_rs::from_reader;


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct RunInfo {
    #[serde(rename = "Version", default)]
    pub version: u32,
    #[serde(rename = "Run")]
    pub runs: Run
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Run {
    #[serde(rename = "Id", default)]
    pub id: String, // parsed from run info, outputed in fastq
    #[serde(rename = "Number", default)]
    pub number: u64,
    #[serde(rename = "Flowcell", default)]
    pub flowcell: String,
    #[serde(rename = "Instrument", default)]
    pub instrument: String,
    #[serde(rename = "Date", default)]
    pub date: String,
    #[serde(rename = "Reads")]
    pub reads: Reads,
    #[serde(rename = "FlowcellLayout")]
    pub flow_cell_layout: FlowcellLayout,
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Reads {
    #[serde(rename = "Read", default)]
    pub read: Vec<Read>,
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Read {
    #[serde(rename = "Number", default)]
    pub number: u64, // parsed from run info, outputed in fastq
    #[serde(rename = "NumCycles", default)]
    pub num_cycles: u64, // number of cycles expected for one read (~100)
    #[serde(rename = "IsIndexedRead", deserialize_with = "bool_from_string")]
    pub is_indexed_read: bool,
}

fn bool_from_string<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: de::Deserializer<'de>,
{
    match String::deserialize(deserializer)?.as_ref() {
        "Y" => Ok(true),
        "N" => Ok(false),
        other => Err(de::Error::invalid_value(
            de::Unexpected::Str(other),
            &"Y or N",
        )),
    }
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct FlowcellLayout {
    #[serde(rename = "LaneCount", default)]
    pub lane_count: u32,
    #[serde(rename = "SurfaceCount", default)]
    pub surface_count: u32,
    #[serde(rename = "SwathCount", default)]
    pub swath_count: u64,
    #[serde(rename = "TileCount", default)]
    pub tile_count: u64,
    #[serde(rename = "FlowcellSide", default)]
    pub flowcell_side: u32,
    #[serde(rename = "TileSet")]
    pub tile_set: TileSet,
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct TileSet {
    #[serde(rename = "TileNamingConvention", default)]
    pub tile_naming_convention: String,
    #[serde(rename = "Tiles")]
    pub tiles: Tiles
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Tiles {
     #[serde(rename = "Tile", default)]
    pub tile: Vec<String>,
}


/// Parse a `RunInfo.xml` file into a `RunInfo` struct or panic
pub fn parse_run_info(run_info_path: &Path) -> RunInfo {
    let run_xml = match fs::read_to_string(run_info_path) {
        Err(e) => panic!(
            "couldn't read {}: {}", run_info_path.display(), e
        ),
        Ok(s) => s,
    };

    let runinfo: RunInfo = match from_reader(run_xml.as_bytes()) {
        Err(e) => panic!(
            "Error parsing RunInfo: {}", e
        ),
        Ok(r) => r,
    };

    runinfo
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_runinfo() {
        let filename_info = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/RunInfo.xml");
        let actual_runinfo = parse_run_info(filename_info);
        let expected_runinfo =
            RunInfo {
                version: 5,
                runs: Run {
                    id: "190414_A00111_0296_AHJCWWDSXX".to_owned(),
                    number: 296,
                    flowcell: "HJCWWDSXX".to_owned(),
                    instrument: "A00111".to_owned(),
                    date: "4/14/2019 1:17:20 PM".to_owned(),
                    reads: Reads {
                        read: vec![
                            Read { number: 1, num_cycles: 4, is_indexed_read: false },
                            Read { number: 2, num_cycles: 8, is_indexed_read: true },
                            Read { number: 3, num_cycles: 8, is_indexed_read: true },
                            Read { number: 4, num_cycles: 4, is_indexed_read: false },
                        ]
                    },
                    flow_cell_layout: FlowcellLayout {
                        lane_count: 1,
                        surface_count: 1,
                        swath_count: 6,
                        tile_count: 3,
                        flowcell_side: 1,
                        tile_set: TileSet {
                            tile_naming_convention: "FourDigit".to_owned(),
                            tiles: Tiles { 
                                tile: vec![
                                    "1_1101".to_owned(),
                                    "1_1102".to_owned(),
                                    "1_1103".to_owned(),
                                ]
                            }
                        },
                    }
                }
            };
        assert_eq!(actual_runinfo, expected_runinfo)
    }

    #[test]
    #[should_panic(
      expected = r#"couldn't read test_data/no_RunInfo.xml: No such file or directory (os error 2)"#
    )]
    fn test_runinfo_no_file() {
        let filename_info = Path::new("test_data/no_RunInfo.xml");
        parse_run_info(filename_info);
    }

    #[test]
    #[should_panic(
      expected = r#"Error parsing RunInfo: 4:10 Unexpected closing tag: RunInfo, expected Run"#
    )]
    fn test_runinfo_bad_file() {
        let filename_info = Path::new("test_data/bad_RunInfo.xml");
        parse_run_info(filename_info);
    }
}
