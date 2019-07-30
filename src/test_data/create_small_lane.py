import pathlib
import struct
import numpy as np
import os

def create_small_lane(test_path, test_dir_name, lane_path, wanted_num_tiles):
    # initialize subpath storage
    filter_paths = []
    cbcl_paths = []

    # open lane_path and separate files (filters) from directories (cbcl_paths)
    os.chdir(pathlib.Path(test_path))
    for subpath in os.listdir(pathlib.Path(lane_path)):
        if os.path.isfile(os.path.join(lane_path, subpath)):
            if len(filter_paths < wanted_num_tiles):
                filter_paths.append(subpath)
            else:
                break
        else:
            cbcl_paths.append(subpath)

    # create test_lane directory
    small_lane_path = pathlib.Path(os.path.join(test_path, test_dir_name))
    os.mkdir(small_lane_path)
    os.chdir(small_lane_path)


    # iterate through cbcl_paths and for each one:
        # 1) make dir for the new cbcl path
        # 2) iterate through the tile_offsets to determine number of lane parts needed
        # 3) make dir for however many lane parts are necessary
        # 4) copy the bytes for length of cbcl header + one tile for the respective lane parts

    for cbcl_path in cbcl_paths:
        os.chdir(small_lane_path)
        small_cbcl_path = os.path.join(small_lane_path, cbcl_path)
        os.mkdir(small_cbcl_path)
        os.chdir(small_cbcl_path)
        full_cbcl_path = os.path.join(lane_path, cbcl_path)

        for lane_part_path in os.listdir(full_cbcl_path):  # only need to go through first lane_part if making a small file (break after first loop)
            small_lane_part_path = os.path.join(small_cbcl_path, lane_part_path)
            full_lane_part_path = pathlib.Path(os.path.join(full_cbcl_path, lane_part_path))

            with full_lane_part_path.open(mode="rb") as f:
                print(full_lane_part_path)
                first = f.read(2)  # 1
                header_size = struct.unpack("<I", f.read(4))
                packed_header_size = struct.pack("<I", 49 + wanted_num_tiles*16)  # 2
                third = f.read(38)  # skip to number of tiles, 3
                packed_num_tiles = struct.pack("<I", wanted_num_tiles)  # 4
                actual_num_tiles = struct.unpack("<I", f.read(4))[0]
                tile_offsets_bytes = f.read(wanted_num_tiles * 16)  # 5
                _ = f.read((actual_num_tiles - wanted_num_tiles) * 16)  # skip the rest of the tile_offsets
                non_PF_clusters_excluded = f.read(1)  # 6

                # read the first "wanted_num_tiles" number of tile_offsets
                tile_offsets = []
                for tile_num in range(wanted_num_tiles):
                    # check the shape of this thing you're appending!!
                    tile_offsets.append(struct.unpack("<IIII", f.read(16)))
                print(tile_offsets)

                tiles = f.read(sum([offset[3] for offset in tile_offsets]))  # read in the sum of all the compressed block sizes of desired tiles, 7

            os.mkdir(small_lane_part_path)
            os.chdir(small_lane_part_path)
            small_lane_part = open("small_lane_part.cbcl", "wb")
            small_lane_part.write(first)  # 1
            small_lane_part.write(packed_header_size)  # 2
            small_lane_part.write(third)  # 3
            small_lane_part.write(packed_num_tiles)  # 4
            small_lane_part.write(tile_offsets_bytes)  # 5
            small_lane_part.write(non_PF_clusters_excluded)  # 6
            small_lane_part.write(tiles)  # 7
            small_lane_part.close()

            break

# create_small_lane("C:/Users/gabby.shvartsman/code/bcl2fastr/src/test_data", "test_lane", pathlib.Path("C:/Users/gabby.shvartsman/code/bcl2fastr/src/test_data/L001"), 1)
