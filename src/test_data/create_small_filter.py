import pathlib
import struct

def create_small_filter(filter_file, wanted_num_clusters):
    with filter_file.open(mode = "rb") as f:
        first = f.read(8) # first part of the file: _, version
        packed_num_clusters = struct.pack("<I", wanted_num_clusters) # second part of the file: num_clusters
        _ = f.read(4) # skip the actual num_clusters
        bools = f.read(wanted_num_clusters) # third part of the file: bool cluster values

    test_filter = open("test_filter.filter", "wb")
    test_filter.write(first)
    test_filter.write(packed_num_clusters)
    test_filter.write(bools)
    test_filter.close()


if __name__ == "__main__":
    create_small_filter(pathlib.Path("s_1_1101.filter"), 5)
