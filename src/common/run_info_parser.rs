//! Deserializes the `RunInfo.xml` file from a NovaSeq run into a useful struct
//! of information about the sequencing run.

use std::{
    fs::File,
    path::Path,
};
use serde::{de, Deserialize};
use serde_xml_rs::from_reader;


/// The top-level struct for the contents of RunInfo.xml
#[derive(Debug, PartialEq, Eq)]
pub struct RunInfo {
    /// Version number of this file (depends on the sequencer)
    pub version: u32,
    /// Full instrument id string (date, instrument, number, flowcell)
    pub id: String,
    /// Number representing how many runs this instrument has performed
    pub number: u64,
    /// Flowcell serial number
    pub flowcell: String,
    /// Instrument serial number/identifier
    pub instrument: String,
    /// The date and time of the run
    pub date: String,
    /// Format of the run: number of reads, read lengths, and which are indexes
    pub reads: Vec<Read>,
    /// Flowcell information: number of lanes, surfaces, and tiles
    pub flowcell_layout: FlowcellLayout,
}


/// Deserialize RunInfo, including flattening the inner Run struct
/// into the top level
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

        #[derive(Deserialize)]
        struct Reads {
            #[serde(rename = "Read")]
            pub read: Vec<Read>,
        }

        fn reads_to_vec<'de, D>(deserializer: D) -> Result<Vec<Read>, D::Error>
        where
            D: de::Deserializer<'de>,
        {
            let reads = Reads::deserialize(deserializer)?;

            let reads = reads.read.into_iter().scan(0, |i, r| {
                *i += r.num_cycles;
                Some(Read { start: *i - r.num_cycles, end: *i, .. r })
            }).collect();

            Ok(reads)
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


/// Information about one of the reads in a run
#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Read {
    /// Which read this is
    #[serde(rename = "Number")]
    pub number: u64,
    /// How many cycles (e.g. bases) in the read
    #[serde(rename = "NumCycles")]
    pub num_cycles: usize,
    /// Whether or not it is an index read
    #[serde(rename = "IsIndexedRead", deserialize_with = "bool_from_string")]
    pub is_indexed_read: bool,
    /// The starting index of the read within the full base array.
    /// This is not in the XML file but is calculated during deserialization
    #[serde(default)]
    pub start: usize,
    /// The end index of the read within the full base array.
    /// This is not in the XML file but is calculated during deserialization
    #[serde(default)]
    pub end: usize,
}


/// Convert from Y or N character to a boolean
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


/// Information about the flowcell used in the run
#[derive(Debug, PartialEq, Eq)]
pub struct FlowcellLayout {
    /// Number of lanes
    pub lane_count: usize,
    /// Number of surfaces per lane
    pub surface_count: usize,
    /// Swathes per lane
    pub swath_count: u64,
    /// Number of tiles per swath
    pub tile_count: u64,
    /// Sides of the flowcell
    pub flowcell_side: u32,
    /// Format for naming tiles
    pub tile_naming_convention: String,
    /// A Vec of tile names
    pub tiles: Vec<String>,
}


/// Deserialize the FlowcellLayout struct including flattening the interior
/// TileNamingConvention struct into the top level
impl<'de> Deserialize<'de> for FlowcellLayout {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>
    {
        #[derive(Deserialize)]
        struct Outer {
            #[serde(rename = "LaneCount")]
            lane_count: usize,
            #[serde(rename = "SurfaceCount")]
            surface_count: usize,
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
            let tiles = Tiles::deserialize(deserializer)?;

            Ok(tiles.tile)
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
    let run_xml = File::open(run_info_path)?;

    let runinfo: RunInfo = match from_reader(run_xml) {
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
                    Read { number: 1, start: 0, end: 4, num_cycles: 4, is_indexed_read: false },
                    Read { number: 2, start: 4, end: 12, num_cycles: 8, is_indexed_read: true },
                    Read { number: 3, start: 12, end: 20, num_cycles: 8, is_indexed_read: true },
                    Read { number: 4, start: 20, end: 24, num_cycles: 4, is_indexed_read: false },
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
    fn no_reads() {
        let filename_info = Path::new("test_data/bad_RunInfo_no_reads.xml");
        parse_run_info(filename_info).unwrap();
    }

    #[test]
    #[should_panic(
        expected = r#"Error parsing RunInfo: custom: 'missing field `Tile`'"#
    )]
    fn no_tiles() {
        let filename_info = Path::new("test_data/bad_RunInfo_no_tiles.xml");
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
