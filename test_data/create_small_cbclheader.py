import pathlib
import struct
import numpy as np

def create_small_cbclheader(cbcl_file, wanted_num_tiles):
    with cbcl_file.open(mode="rb") as f:
        first = f.read(44) # first part of the file: version, header_size, bits_, bits_, num_bins, bins
        number_of_tile_records = wanted_num_tiles
        packed_num_tiles = struct.pack("<I", number_of_tile_records) # second part of file: number_of_tile_records
        actual_num_tiles = struct.unpack("<I", f.read(4))[0] # read the actual number_of_tile_records
        tile_offsets = f.read(wanted_num_tiles * 16) # third part of the file: first 10 tile_offsets (4x4x10 = 160 (each part is I))
        _ = f.read((actual_num_tiles - wanted_num_tiles) * 16) # skip the rest of the tile_offsets

        non_PF_clusters_excluded = f.read(1) # last part of the file: non_PF_clusters_excluded
        check = struct.unpack("<B", non_PF_clusters_excluded)[0]

    test_cbcl = open("test_cbcl_header.cbcl", "wb")
    test_cbcl.write(first)
    test_cbcl.write(packed_num_tiles)
    test_cbcl.write(tile_offsets)
    test_cbcl.write(non_PF_clusters_excluded)
    test_cbcl.close()


if __name__ == "__main__":
    create_small_cbclheader(pathlib.Path("L001_1.cbcl"), 5)
