//defining structs
struct RunInfo {
    run_path: String, // path to specific run dir
    run_id: String, // parsed from run info, outputed in fastq
    num_cycles: u32, // number of cycles
    num_lanes: u32, // number of lanes on the cell
    len_read: u32, // in bp (usally ~100)
    len_index: u32, // in bp (usually 8 or 16)
    tile_coords: Vec<(u32, u32)>, // (x, y) coords of tiles in each cbcl file (same in each file)
    tile_dims: (u32, u32), // how many tiles on x vs y axis for each cbcl file (same in each file)
    cbcls: Vec<CBCL>, // number of cbcl files depends on num_cycles
    qscore_bins: Vec<(f64, f64)>, // (min, max) of range for corresponding bins to simplify raw qscores to
    read_n: u32, // number of reads to expect
}

struct CBCL {
    lane: i32, // which lane the file is from
    lane_part: i32, // which lane_part file it's from
    num_tiles: i32, // num_tiles = RunInfo.title_dims.0*RunInfo.title_dims.1
    tiles: [Tile; num_tiles],
}

struct Tile {
    bases: [[i32; RunInfo.tile_dims.1]; RunInfo.tile_dims.0], // matrix of bases per tile; want to convert to the 0, 1, 2, 3, 4 notation in the process of building a Tile struct
    qscores: [[i32; RunInfo.tile_dims.1]; RunInfo.tile_dims.0], // matrix of qscores; want to do the raw qscore --> bin conversion when building a Tile struct
    filter: , //dims?
}

// reading in/building structs
fn parse_run_info() {}
fn parse_cbcl() {}
fn parse_tile() {}


// structs tests
#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_runinfo {
        // do we need to check types here? because won't the compiler already error if you're building a struct with an incorrectly-typed attribute?
    }

}

// INCLUDE IN TESTING:
// make sure to include if cbcl files = None --> would be output by sequencer if clustering is super bad and the files were just super messed up
// if the run is complete, there should be the same number of bcl files for each lane so check for that
// match the number of cbcl files expected based on the num_cycles specified (check with Michelle and Rene on the exact file structure)
