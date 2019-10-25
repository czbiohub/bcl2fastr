#!/usr/bin/env python

import gzip
import itertools
import pathlib
import struct
import xml.etree.ElementTree as et

from typing import Optional, Mapping

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np


@dataclass(eq=False, frozen=True)
class CBCLHeader:
    # probably all constant, but tiny
    version: int
    size: int
    bits_per_basecall: int
    bits_per_qscore: int
    number_of_bins: int
    bins: np.ndarray
    number_of_tile_records: int
    # variable
    tile_offsets: np.ndarray
    non_PF_clusters_excluded: bool


get_cycle = lambda cfn: int(cfn.parent.name[1:-2])
get_part = lambda cfn: int(cfn.name[2])
get_tile = lambda cfn: int(cfn.name[4:8])


class NovaSeqRun(object):
    LANES = range(1, 5)  # supports 4 lanes per run
    PARTS = range(1, 3)  # supports 2 parts per lane

    def __init__(
        self,
        run_path: pathlib.Path,
        n_threads: int = 8,
        start: int = 1,
        stop: Optional[int] = None,
    ):
        # path to specific run directory
        self.run_path = run_path

        # xy coordinates for the location of each read per tile
        self.locs = self.read_locs(run_path / "Data" / "Intensities" / "s.locs")

        # generate the naming format used in fastq files
        _, instrument, run_no, flowcell_id = run_path.name.split("_")
        self.run_id = f"@{instrument}:{int(run_no)}:{flowcell_id}"

        self.run_info = self.read_run_info(run_path)
        self.read_n = range(sum(1 for i1, i2, index in self.run_info if not index))

        self.cbcl_files = dict()
        self.filters = dict()
        self.headers = dict()
        self.tiles = dict()

        # assign each tile to a separate thread
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            for lane, part, cbcl_files in executor.map(
                NovaSeqRun.get_file_lists,
                itertools.repeat(run_path),
                *zip(*itertools.product(self.LANES, self.PARTS)),
                itertools.repeat(start),
                itertools.repeat(stop),
            ):
                if cbcl_files:
                    self.cbcl_files[lane, part] = cbcl_files

            self.lane_parts = sorted(self.cbcl_files)

            for lane, filters in executor.map(
                NovaSeqRun.read_filters, itertools.repeat(run_path), self.LANES
            ):
                if filters:
                    self.filters[lane] = filters

            for lane, part, headers, tiles in executor.map(
                NovaSeqRun.read_headers,
                *zip(
                    *(
                        (lane, part, self.cbcl_files[lane, part])
                        for lane, part in self.lane_parts
                    )
                ),
            ):
                self.headers[lane, part] = headers
                self.tiles[lane, part] = tiles

        self.code_set = self.get_code_set(self.headers)

        cycles = {range(min(v), max(v) + 1) for v in self.cbcl_files.values()}
        assert len(cycles) == 1

        self.cycles = cycles.pop()

    # look at a subset of every read (ex: you want only the indices or only the reads)
    def subset_cycles(self, start: int, stop: int):
        if start < 1:
            raise ValueError("start must be ≥ 1")
        if stop > len(self.cbcl_files[self.lane_parts[0]]) + 1:
            raise ValueError("'stop' must be ≤ the number of cycles")
        self.cycles = range(start, stop)

    # use run_info to get only the read portions of a sequence
    def get_reads(self, seq):
        return tuple(seq[i1:i2] for i1, i2, index in self.run_info if not index)

    # use run_info to get only the index portions of a sequence
    def get_indexes(self, seq):
        return tuple(seq[i1:i2] for i1, i2, index in self.run_info if index)

    # parse through the run meta-info
    @staticmethod
    def read_run_info(run_path: pathlib.Path):
        with (run_path / "RunInfo.xml").open() as f:
            run_info = et.parse(f)

        read_elems = run_info.findall("./Run/Reads/Read[@NumCycles][@Number]")
        read_elems.sort(key=lambda el: int(el.get("Number")))

        read_info = [(0, 0, False)]

        for i, elem in enumerate(read_elems):
            last_i = read_info[i][1]
            read_info.append(
                (
                    last_i,
                    last_i + int(elem.get("NumCycles")),
                    elem.get("IsIndexedRead") == "Y",
                )
            )

        return read_info[1:]

    # get_file_lists: store list of CBCL files for a given lane + part
    @staticmethod
    def get_file_lists(
        run_path: pathlib.Path,
        lane: int,
        part: int,
        start: int,
        stop: Optional[int] = None,
    ):
        if stop is None:
            stop = np.inf

        cbcl_file_list = [
            cbcl_file
            for cbcl_file in (
                run_path / "Data" / "Intensities" / "BaseCalls" / f"L00{lane}"
            ).glob(f"C*.1/L00{lane}_{part}.cbcl")
            if start <= get_cycle(cbcl_file) < stop
        ]
        if cbcl_file_list:
            return lane, part, {get_cycle(cfn): cfn for cfn in cbcl_file_list}
        else:
            return lane, part, dict()

    # read_filters: read in the filter data for each tile of a given lane
    @staticmethod
    def read_filters(run_path: pathlib.Path, lane: int):
        filter_list = list(
            (run_path / "Data" / "Intensities" / "BaseCalls" / f"L00{lane}").glob(
                f"s_{lane}_*.filter"
            )
        )

        if filter_list:
            return (
                lane,
                {get_tile(ffn): NovaSeqRun.read_filter(ffn) for ffn in filter_list},
            )
        else:
            return lane, dict()

    # read_headers: read in the CBCL header information for a given lane + part
    @staticmethod
    def read_headers(lane: int, part: int, cbcl_files: Mapping[int, pathlib.Path]):
        header_data = {c: NovaSeqRun.read_header(cbcl_files[c]) for c in cbcl_files}

        # assertion: the tile numbers are all the same in this sequence of files
        assert np.array_equal(
            np.vstack([ch.tile_offsets[:, 0] for ch in header_data.values()]).max(0),
            np.vstack([ch.tile_offsets[:, 0] for ch in header_data.values()]).min(0),
        )

        tiles = next(iter(header_data.values())).tile_offsets[:, 0]

        return lane, part, header_data, tiles

    # read_filter: load a single array of filter values
    @staticmethod
    def read_filter(filter_file: pathlib.Path):
        with filter_file.open(mode="rb") as f:
            _, filter_version, n_clusters = struct.unpack("<III", f.read(12))
            pf = np.fromfile(f, dtype=np.uint8, count=n_clusters)

        return (pf & 0b1).astype(bool)

    # save new CBCLHeader object to easily access header info from input file
    @staticmethod
    def read_header(cbcl_file: pathlib.Path):
        with cbcl_file.open(mode="rb") as f:
            version, header_size, bits_per_basecall, bits_per_qscore, number_of_bins = struct.unpack(
                "<HIBBI", f.read(12)
            )
            # bins compress a wide range of qscores into 4 scores (0, 1, 2, 3)
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

            cbcl_header = CBCLHeader(
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

        return cbcl_header

    # takes set of characters for the ranges of bins and decodes as ASCII characters
    # which is the proper output format for fastq
    @staticmethod
    def get_code_set(headers):
        # retrieve qscore binning
        code_set = set(
            tuple(c_h.bins[:, 1] + 33)
            for lane_part in headers
            for c_h in headers[lane_part].values()
        )

        # assertion: bins are constant across a whole run
        assert len(code_set) == 1
        return bytes(code_set.pop()).decode() + "#"

    # read the location file for the run (fixed for all lanes)
    @staticmethod
    def read_locs(loc_file: pathlib.Path):
        with loc_file.open("rb") as f:
            _, _, n_clusters = struct.unpack("<III", f.read(12))
            cbcl_locs = np.fromfile(f, dtype=np.float32, count=2 * n_clusters)

        # convert to the float type that matches those in fastq
        cbcl_locs = (cbcl_locs.astype(float) * 10 + 1000).reshape((n_clusters, 2))

        return cbcl_locs.round().astype(int)


# parses data based on the byte-by-byte breakdown in bcl2fastq2 software guide
def get_byte_lists(
    novaseq_run: NovaSeqRun, lane: int, part: int, i: int, cf: np.ndarray
):  # cf is the cluster-pass filter for a specific tile, i is the tile number

    # note if an odd number of reads pass the filter since then when parsing by
    # bytes, the last half-byte does not represent any real information
    odd = -(cf.sum() % 2) or None

    for c in novaseq_run.cycles:
        c_h = novaseq_run.headers[lane, part][c]

        with novaseq_run.cbcl_files[lane, part][c].open("rb") as f:
            f.seek(c_h.size + c_h.tile_offsets[:i, 3].sum(dtype=int))

            try:
                byte_array = np.frombuffer(
                    gzip.decompress(f.read(c_h.tile_offsets[i, 3])),
                    dtype=np.uint8,
                    count=c_h.tile_offsets[i, 2],
                )
            except OSError:
                yield (None, None)
                continue

            unpacked = np.unpackbits(byte_array[::-1]).reshape((-1, 2, 2))
            repacked = (np.packbits(unpacked, axis=2) >> 6).squeeze()[::-1]

            # a set of tuples for each lane part that includes (base_array, qscore_array)
            # for each tile in that lane part
            if c_h.non_PF_clusters_excluded:
                yield repacked[:odd, 1], repacked[:odd, 0]
            else:
                yield repacked[cf, 1], repacked[cf, 0]


# extract the reads of a tile in a specific lane and lane part
def extract_reads(novaseq_run: NovaSeqRun, lane: int, part: int, i: int):
    tile_no = novaseq_run.tiles[lane, part][i]
    cf = novaseq_run.filters[lane][tile_no]

    # set matrices of bases and qscores across cycles
    base_matrix = 4 * np.ones((cf.sum(), len(novaseq_run.cycles)), dtype=np.uint8)
    qscore_matrix = 4 * np.ones_like(base_matrix)

    ba_generator = enumerate(get_byte_lists(novaseq_run, lane, part, i, cf))

    # bases: 0, 1, 2, or 3 (corresponding to each base)
    # qscores: 0, 1, 2, or 3 (corresponding to each qscore bin)
    for j, (base_array, qscore_array) in ba_generator:
        if base_array is not None:
            base_matrix[:, j] = base_array
            qscore_matrix[:, j] = qscore_array

    # if quality score is 0 it gets masked to the minimum value
    base_matrix[qscore_matrix == 0] = 4
    qscore_matrix[qscore_matrix == 0] = 4

    # return a tuple of (bases, qscores, title of tile) for each read in a tile
    yield from zip(
        (
            "".join("ACGTN"[b] for b in base_matrix[k, :])
            for k in range(base_matrix.shape[0])
        ),
        (
            "".join(novaseq_run.code_set[b] for b in qscore_matrix[k, :])
            for k in range(qscore_matrix.shape[0])
        ),
        (
            f"{novaseq_run.run_id}:{lane}:{tile_no}:{r[0]}:{r[1]}"
            for r in novaseq_run.locs[cf, :]
        ),
    )
