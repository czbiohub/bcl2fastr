

//defining structs
struct RunInfo {
    run_path: String, // path to specific run dir
    run_id: String, // parsed from run info, outputed in fastq
    paired_end: Boolean, // if paired-end read --> True, else False (determines if you double cycle numbers)
    num_read_cycles: u32, // number of cycles expected for one read (~100)
    num_index_cycles: u32, // number of cycles expected for one index (~8)
    num_lanes: u32, // number of lanes on the cell
    tile_coords: Vec<(u32, u32)>, // (x, y) coords of tiles in each cbcl file (same in each file)
    tile_dims: (u32, u32), // how many tiles on x vs y axis for each cbcl file (same in each file)
    qscore_bins: Vec<(f64, f64)>, // (min, max) of range for corresponding bins to simplify raw qscores to
    read_n: u32, // number of reads to expect
    lanes: Vec<Lane>, // store each lane (will parse by lane)
}

struct Lane {
    cbcls: Vec<CBCL>, // cbcl files for that particular lane
    sample: u32, // which sample this lane is a part of (to make sure to differentiate in xp runs)
}

struct CBCL {
    lane: u32, // which lane the file is from
    lane_part: u32, // which lane_part file it's from
    num_tiles: u32, // num_tiles = RunInfo.title_dims.0*RunInfo.title_dims.1
    tiles: [Tile; num_tiles],
}

struct Tile {
    bases: [[u32; RunInfo.tile_dims.1]; RunInfo.tile_dims.0], // matrix of bases per tile; want to convert to the 0, 1, 2, 3, 4 notation in the process of building a Tile struct
    qscores: [[u32; RunInfo.tile_dims.1]; RunInfo.tile_dims.0], // matrix of qscores; want to do the raw qscore --> bin conversion when building a Tile struct
    filter: , //dims?
}

// reading in/building a RunInfo object given the runpath
fn parse_run_info(run_path: String) -> RunInfo {
    //
    // return RunInfo {
    //     run_path: run_path,
    //     run_id: ,
    //     paired_end: ,
    //     num_read_cycles: ,
    //     num_index_cycles: ,
    //     num_lanes: ,
    //     tile_coords: ,
    //     tile_dims: ,
    //     qscore_bins: ,
    //     read_n: ,
    //     lanes: ,
    // }

}

fn parse_lane(lane_path: String) -> Lane {

}

fn parse_cbcl(cbcl_path: String) -> CBCL {

}

fn parse_tile(tile_path: String) -> Tile {

}


// structs tests
#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    // each lane should have the same number of cbcl files as the expected number of cycles in the run
    fn lanes_match_cycles(run_path: String) {
        file_path: String = run_path; // insert file path of corrupted file
        run: RunInfo = parse_run_info(file_path);
        expected_cbcls = 0;
        if run.paired_end {
            expected_cbcls = (run.num_read_cycles + run.num_index_cycles) * 2;
        } else {
            expected_cbcls = run.num_read_cycles + run.num_index_cycles;
        }
        for lane_num in 1..run.lanes.len() {
            actual_cbcls = run.lanes[lane_num].cbcls.len()
            assert_eq!(expected_cbcls, actual_cbcls, "Expected Lane {} to have {} CBCL files but it had {}", lane_num, expected_cbcls, actual_cbcls);
        }
    }

}
