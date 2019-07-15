# import pathlib
import struct
#
# from dataclasses import dataclass

import numpy as np

#
# @dataclass(eq=False, frozen=True)
# class CBCLHeader:
#     # probably all constant, but tiny
#     version: int
#     size: int
#     bits_per_basecall: int
#     bits_per_qscore: int
#     number_of_bins: int
#     bins: np.ndarray
#     number_of_tile_records: int
#     # variable
#     tile_offsets: np.ndarray
#     non_PF_clusters_excluded: bool
#

def read_header(cbcl_file):
    with cbcl_file.open(mode="rb") as f:
        version, header_size, bits_per_basecall, bits_per_qscore, number_of_bins = struct.unpack(
            "<HIBBI", f.read(12)
        )
        #bins compress a wide range of qscores into 4 scores (0, 1, 2, 3)
        bins = (
            np.fromfile(f, dtype=np.uint32, count=2 * number_of_bins)
            .reshape((number_of_bins, 2))
            .astype(np.uint8)
        )

        number_of_tile_records = struct.unpack("<I", f.read(4))[0]
        tile_offsets = np.fromfile(
            f, dtype=np.uint32, count=4 * number_of_tile_records
        ).reshape((number_of_tile_records, 4))

        non_PF_clusters_excluded = bool(struct.unpack("B", f.read(1))[0])

        println(
            version,
            header_size,
            bits_per_basecall,
            bits_per_qscore,
            number_of_bins,
            bins,
            number_of_tile_records,
            tile_offsets,
            non_PF_clusters_excluded,
        )

    # return cbcl_header

def main():
    read_header("L001_1_cbcl_header.cbcl")
