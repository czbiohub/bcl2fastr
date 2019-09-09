import pathlib
import struct

def create_small_locs(locs_file, wanted_num_clusters):
    with locs_file.open(mode="rb") as f:
        first = f.read(8) # first part of the file: _ and num_clusters
        num_clusters = wanted_num_clusters
        packed_num_clusters = struct.pack("<I", num_clusters) # second part of the file: num_clusters
        _ = f.read(4) # skip the actual num_clusters in the file
        locs = f.read(wanted_num_clusters * 8) # third part of the file: only "wanted_num_clusters" number of locs

    test_locs = open("test_locs.locs", "wb")
    test_locs.write(first)
    test_locs.write(packed_num_clusters)
    test_locs.write(locs)
    test_locs.close()


if __name__ == "__main__":
    create_small_locs(pathlib.Path("s.locs"), 5)
