import pathlib
import struct
import os
import glob
import shutil
import sys


def create_small_lane(test_path, lane_path, wanted_num_tiles):
    # initialize subpath storage
    # iterate through cbcl_paths and for each one:
        # 1) make dir for the new cbcl path
        # 2) iterate through the tile_offsets to determine number of lane parts needed
        # 3) make dir for however many lane parts are necessary
        # 4) copy the bytes for length of cbcl header + one tile for the respective lane parts
    
    test_lane_path = os.path.join(test_path, "test_lane")

    if os.path.exists(test_lane_path):
        shutil.rmtree(test_lane_path)
        
    os.mkdir(test_lane_path)
    
    filter_paths = sorted(
        glob.glob(os.path.join(
            lane_path,
            "*filter"))
    )[0:wanted_num_tiles]
        
    for filter_path in filter_paths:
        shutil.copyfile(filter_path, os.path.join(test_lane_path, os.path.basename(filter_path)))
 
    lane_name = os.path.basename(lane_path)
    actual_cbcl_files = sorted(glob.glob(os.path.join(lane_path, "C*", f"{lane_name}_1.cbcl")))
 
    for cbcl_file in actual_cbcl_files:
        # test_dir_cbcl_path:  /usr/src/bcl2fastr/src/test_data/test_lane/C1.1
        # test_lane_basename:  /usr/src/bcl2fastr/src/test_data/test_lane/C1.1/L001_1.cbcl
        test_dir_cbcl_path = os.path.join(test_lane_path, os.path.basename(os.path.dirname(cbcl_file)))
        test_lane_basename = os.path.join(test_dir_cbcl_path, os.path.basename(cbcl_file))
        
        os.mkdir(test_dir_cbcl_path)
        
        with open(cbcl_file, mode="rb") as f:
            first = f.read(2)  # 1
            _ = f.read(4)  # skipping the header size
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
                tile_offsets.append(
                    struct.unpack("<IIII", tile_offsets_bytes[tile_num * 16 : (tile_num + 1) * 16])
                )
                
            print(tile_offsets)

            # read in the sum of all the compressed block sizes of desired tiles, 7
            tiles = f.read(sum([offset[3] for offset in tile_offsets]))  

        with open(test_lane_basename, "wb") as f:
            f.write(first)  # 1
            f.write(packed_header_size)  # 2
            f.write(third)  # 3
            f.write(packed_num_tiles)  # 4
            f.write(tile_offsets_bytes)  # 5
            f.write(non_PF_clusters_excluded)  # 6
            f.write(tiles)  # 7

            
if __name__ == "__main__":
    save_path = sys.argv[1]
    data_path = sys.argv[2]
    create_small_lane(save_path, data_path, 1)

