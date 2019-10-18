//! Deserializes the `RunInfo.xml` file from a NovaSeq run into a useful struct
//! of information about the sequencing run.

use std::{
    fs,
    path::Path,
};
use serde::{de, Deserialize};
use serde_xml_rs::from_reader;


#[derive(Debug, PartialEq, Eq)]
pub struct RunInfo {
    pub version: u32,
    pub id: String,
    pub number: u64,
    pub flowcell: String,
    pub instrument: String,
    pub date: String,
    pub reads: Vec<Read>,
    pub flowcell_layout: FlowcellLayout,
}


impl<'de> Deserialize<'de> for RunInfo {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>
    {
        #[derive(Deserialize)]
        struct Outer {
            #[serde(rename = "Version")]
            version: u32,
            #[serde(rename = "Run")]
            run: Inner
        }

        #[derive(Deserialize)]
        struct Inner {
            #[serde(rename = "Id")]
            id: String,
            #[serde(rename = "Number")]
            number: u64,
            #[serde(rename = "Flowcell")]
            flowcell: String,
            #[serde(rename = "Instrument")]
            instrument: String,
            #[serde(rename = "Date")]
            date: String,
            #[serde(rename = "Reads", deserialize_with = "reads_to_vec")]
            reads: Vec<Read>,
            #[serde(rename = "FlowcellLayout")]
            flowcell_layout: FlowcellLayout,
        }

        let helper = Outer::deserialize(deserializer)?;

        Ok(RunInfo {
            version: helper.version,
            id: helper.run.id,
            number: helper.run.number,
            flowcell: helper.run.flowcell,
            instrument: helper.run.instrument,
            date: helper.run.date,
            reads: helper.run.reads,
            flowcell_layout: helper.run.flowcell_layout,
        })
    }
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Reads {
    #[serde(rename = "Read")]
    pub read: Vec<Read>,
}


fn reads_to_vec<'de, D>(deserializer: D) -> Result<Vec<Read>, D::Error>
where
    D: de::Deserializer<'de>,
{
    let reads = Reads::deserialize(deserializer)?;

    let mut i = 0;
    let reads = reads.read.into_iter().map(|r| {
        let read = Read { start: i, .. r };
        i += r.num_cycles;
        read
    }).collect();

    Ok(reads)
}


#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Read {
    #[serde(rename = "Number")]
    pub number: u64,
    #[serde(rename = "NumCycles")]
    pub num_cycles: usize,
    #[serde(rename = "IsIndexedRead", deserialize_with = "bool_from_string")]
    pub is_indexed_read: bool,
    #[serde(default)]
    pub start: usize,
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


#[derive(Debug, PartialEq, Eq)]
pub struct FlowcellLayout {
    pub lane_count: u32,
    pub surface_count: u32,
    pub swath_count: u64,
    pub tile_count: u64,
    pub flowcell_side: u32,
    pub tile_naming_convention: String,
    pub tiles: Vec<String>,
}


impl<'de> Deserialize<'de> for FlowcellLayout {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>
    {
        #[derive(Deserialize)]
        struct Outer {
            #[serde(rename = "LaneCount")]
            lane_count: u32,
            #[serde(rename = "SurfaceCount")]
            surface_count: u32,
            #[serde(rename = "SwathCount")]
            swath_count: u64,
            #[serde(rename = "TileCount")]
            tile_count: u64,
            #[serde(rename = "FlowcellSide")]
            flowcell_side: u32,
            #[serde(rename = "TileSet")]
            tile_set: Inner,
        }

        #[derive(Deserialize)]
        struct Inner {
            #[serde(rename = "TileNamingConvention")]
            tile_naming_convention: String,
            #[serde(rename = "Tiles", deserialize_with = "tiles_to_vec")]
            tiles: Vec<String>,
        }

        #[derive(Deserialize)]
        struct Tiles {
            #[serde(rename = "Tile")]
            tile: Vec<String>,
        }

        fn tiles_to_vec<'de, D>(deserializer: D) -> Result<Vec<String>, D::Error>
        where
            D: de::Deserializer<'de>,
        {
            match Tiles::deserialize(deserializer) {
                Ok(tiles) => Ok(tiles.tile),
                _ => Err(de::Error::invalid_value(
                    de::Unexpected::Other("Failed to deserialize"),
                    &"a Tiles struct",
                )),
            }
        }

        let helper = Outer::deserialize(deserializer)?;

        Ok(FlowcellLayout {
            lane_count: helper.lane_count,
            surface_count: helper.surface_count,
            swath_count: helper.swath_count,
            tile_count: helper.tile_count,
            flowcell_side: helper.flowcell_side,
            tile_naming_convention: helper.tile_set.tile_naming_convention,
            tiles: helper.tile_set.tiles,
        })
    }
}


/// Parse a `RunInfo.xml` file into a `RunInfo` struct or panic
pub fn parse_run_info(run_info_path: &Path) -> std::io::Result<RunInfo> {
    let run_xml = fs::read_to_string(run_info_path)?;

    let runinfo: RunInfo = match from_reader(run_xml.as_bytes()) {
        Err(e) => panic!(
            "Error parsing RunInfo: {}", e
        ),
        Ok(r) => r,
    };

    Ok(runinfo)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse() {
        let filename_info = Path::new("test_data/190414_A00111_0296_AHJCWWDSXX/RunInfo.xml");
        let actual_runinfo = parse_run_info(filename_info).unwrap();
        let expected_runinfo =
            RunInfo {
                version: 5,
                id: "190414_A00111_0296_AHJCWWDSXX".to_owned(),
                number: 296,
                flowcell: "HJCWWDSXX".to_owned(),
                instrument: "A00111".to_owned(),
                date: "4/14/2019 1:17:20 PM".to_owned(),
                reads: vec![
                    Read { number: 1, start: 0, num_cycles: 4, is_indexed_read: false },
                    Read { number: 2, start: 4, num_cycles: 8, is_indexed_read: true },
                    Read { number: 3, start: 12, num_cycles: 8, is_indexed_read: true },
                    Read { number: 4, start: 20, num_cycles: 4, is_indexed_read: false },
                ],
                flowcell_layout: FlowcellLayout {
                    lane_count: 1,
                    surface_count: 1,
                    swath_count: 6,
                    tile_count: 3,
                    flowcell_side: 1,
                    tile_naming_convention: "FourDigit".to_owned(),
                    tiles: vec![
                        "1_1101".to_owned(),
                        "1_1102".to_owned(),
                        "1_1103".to_owned(),
                    ]
                }
            };
        assert_eq!(actual_runinfo, expected_runinfo)
    }

    #[test]
    #[should_panic(
      expected = r#"No such file or directory"#
    )]
    fn no_file() {
        let filename_info = Path::new("test_data/no_RunInfo.xml");
        parse_run_info(filename_info).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"invalid value: string "Q", expected Y or N"#
    )]
    fn weird_file() {
        let filename_info = Path::new("test_data/weird_RunInfo.xml");
        parse_run_info(filename_info).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"Error parsing RunInfo: custom: 'missing field `Read`'"#
    )]
    fn bad_file() {
        let filename_info = Path::new("test_data/bad_RunInfo.xml");
        parse_run_info(filename_info).unwrap();
    }

    #[test]
    #[should_panic(
      expected = r#"1:1 Unexpected end of stream: no root element found"#
    )]
    fn empty_file() {
        let filename_info = Path::new("test_data/empty_file");
        parse_run_info(filename_info).unwrap();
    }
}
