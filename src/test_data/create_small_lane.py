import pathlib
import struct
import numpy as np
import os

def create_small_lane(test_path_name, lane_path, wanted_num_tiles):
    # initialize subpath storage
    filter_paths = []
    cbcl_paths = []

    # open lane_path and separate files (filters) from directories (cbcl_paths)
    with lane_path.open() as f:
        for subpath in os.listdir(lane_path):
            if os.isfile(os.join(lane_path, subpath)):
                if len(filter_paths < wanted_num_tiles):
                    filter_paths.append(subpath)
                else:
                    break
            else:
                cbcl_paths.append(subpath)

    # create test_lane directory
    small_lane_path = os.mkdir(test_path_name)

    # iterate through cbcl_paths and for each one:
        # 1) make dir for the new cbcl path
        # 2) iterate through the tile_offsets to determine number of lane parts needed
        # 3) make dir for however many lane parts are necessary
        # 4) copy the bytes for length of cbcl header + one tile for the respective lane parts
    for cbcl_path in cbcl_paths:
        os.chdir(small_lane_path)
        small_cbcl_path = os.join(test_path_name, cbcl_path)
        os.mkdir(small_cbcl_path)
        os.chdir(small_cbcl_path)
        # cutoff here ^ for first round of testing
