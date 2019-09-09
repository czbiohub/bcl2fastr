#!/usr/bin/env python

import gzip
import pathlib
import struct

from typing import Mapping

import xml.etree.ElementTree as et

import numpy as np

import click

from seqbot.bcl2fastr import NovaSeqRun, CBCLHeader, get_cycle, get_tile


def make_run_info(
    run_path: pathlib.Path, output_path: pathlib.Path, n_tiles: int, n_bases: int
):
    with (run_path / "RunInfo.xml").open() as f:
        run_info = et.parse(f)

    for e in run_info.findall("./Run/Reads/Read"):
        if e.get("IsIndexedRead") == "N":
            e.set("NumCycles", f"{n_bases // 2}")

    fcl = run_info.find("./Run/FlowcellLayout")
    fcl.attrib.update(
        {"LaneCount": "1", "SurfaceCount": "1", "TileCount": f"{n_tiles}"}
    )

    t = fcl.find("./TileSet/Tiles")
    for i, e in enumerate(sorted(t.findall("./Tile"), key=lambda e: e.text)):
        if i >= n_tiles:
            t.remove(e)

    run_info.write(output_path / "RunInfo.xml", encoding="utf-8", xml_declaration=True)


def subset_locs(run_path: pathlib.Path, output_path: pathlib.Path, n_reads: int):
    # read raw locs
    with (run_path / "Data" / "Intensities" / "s.locs").open("rb") as f:
        a, b, _ = struct.unpack("<III", f.read(12))
        cbcl_locs = f.read(8 * n_reads)

    # write locs to new file
    with (output_path / "Data" / "Intensities" / "s.locs").open("wb") as out:
        out.write(struct.pack("<III", a, b, n_reads))
        out.write(cbcl_locs)


def subset_filters(
    run_path: pathlib.Path, data_path: pathlib.Path, n_tiles: int, n_reads: int
):
    # write filter files
    filter_list = sorted(
        (run_path / "Data" / "Intensities" / "BaseCalls" / f"L001").glob(
            f"s_1_*.filter"
        ),
        key=get_tile,
    )[:n_tiles]

    filter_sums = {}

    for filter_file in filter_list:
        tile = get_tile(filter_file)
        filter_data = NovaSeqRun.read_filter(filter_file)[:n_reads]
        filter_sums[tile] = filter_data.sum()

        with (data_path / filter_file.name).open("wb") as out:
            out.write(struct.pack("<III", 0, 3, n_reads))
            out.write(filter_data.tobytes())

    return filter_sums


def make_header(header: CBCLHeader, new_tile_offsets: np.ndarray):
    header_size = 17 + 4 * header.bins.size + 4 * new_tile_offsets.size

    new_header = [
        struct.pack(
            "<HIBBI",
            header.version,
            header_size,
            header.bits_per_basecall,
            header.bits_per_qscore,
            header.number_of_bins,
        ),
        header.bins.astype(np.uint32).tobytes(),
        struct.pack("<I", new_tile_offsets.shape[0]),
        new_tile_offsets.tobytes(),
        struct.pack("B", header.non_PF_clusters_excluded),
    ]

    return b"".join(new_header)


def subset_cycle(
    cbcl_file: pathlib.Path,
    cbcl_header: CBCLHeader,
    cbcl_filter_sums: Mapping[int, np.ndarray],
    n_tiles: int,
    n_reads: int,
):
    """For a single CBCL file, read the first n_reads from the first n_tiles
    and return the gzip-compressed blocks and tile offsets
    """

    new_tiles = []
    new_tile_offsets = []

    with cbcl_file.open("rb") as f:
        f.seek(cbcl_header.size)

        for i in range(n_tiles):
            tile_i = cbcl_header.tile_offsets[i, 0]
            n_pf = cbcl_filter_sums[tile_i]

            byte_array = np.frombuffer(
                gzip.decompress(f.read(cbcl_header.tile_offsets[i, 3])),
                dtype=np.uint8,
                count=cbcl_header.tile_offsets[i, 2],
            )

            unpacked = np.unpackbits(byte_array).reshape((-1, 2, 2))

            if cbcl_header.non_PF_clusters_excluded:
                byte_array = np.packbits(unpacked[:n_pf, :, :].reshape(4 * n_pf, -1))
            else:
                byte_array = np.packbits(
                    unpacked[:n_reads, :, :].reshape(4 * n_reads, -1)
                )

            new_tiles.append(gzip.compress(byte_array.tobytes()))
            new_tile_offsets.append(
                (
                    tile_i,
                    n_pf if cbcl_header.non_PF_clusters_excluded else n_reads,
                    byte_array.size,
                    len(new_tiles[-1]),
                )
            )

    return b"".join(new_tiles), np.array(new_tile_offsets, dtype=np.uint32)


@click.command()
@click.option(
    "-r",
    "--run_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "-o",
    "--output_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("--n_tiles", required=True, type=int)
@click.option("--n_reads", required=True, type=int)
@click.option("--n_bases", required=True, type=int)
def main(run_path: str, output_path: str, n_tiles: int, n_reads: int, n_bases: int):
    # this is just simpler
    assert n_reads % 2 == 0
    assert n_bases % 2 == 0

    run_path = pathlib.Path(run_path)
    output_path = pathlib.Path(output_path)

    run_info = NovaSeqRun.read_run_info(run_path)

    ix_start = min(i1 for i1, i2, is_index in run_info if is_index)
    ix_end = max(i2 for i1, i2, is_index in run_info if is_index)

    # use n_bases / 2 from start so that we have some unfiltered files
    cycles_to_use = dict(zip(range(1, n_bases // 2 + 1), range(1, n_bases // 2 + 1)))

    # then add the indexes and another n_bases / 2 from second read
    for c in range(ix_start + 1, ix_end + n_bases // 2 + 1):
        cycles_to_use[c] = len(cycles_to_use) + 1

    cbcl_files = {
        get_cycle(cbcl_file): cbcl_file
        for cbcl_file in (
            run_path / "Data" / "Intensities" / "BaseCalls" / "L001"
        ).glob(f"C*.1/L001_1.cbcl")
        if get_cycle(cbcl_file) in cycles_to_use
    }

    _, _, headers, _ = NovaSeqRun.read_headers(1, 1, cbcl_files)

    # write the new RunInfo file
    make_run_info(run_path, output_path, n_tiles, n_bases)

    data_path = output_path / "Data" / "Intensities" / "BaseCalls" / "L001"
    data_path.mkdir(parents=True, exist_ok=True)

    subset_locs(run_path, output_path, n_reads)

    filter_sums = subset_filters(run_path, data_path, n_tiles, n_reads)

    # write cycles
    for c in cycles_to_use:
        (data_path / f"C{cycles_to_use[c]}.1").mkdir(exist_ok=True)

        # take from lane 1 part 1
        new_tiles, new_tile_offsets = subset_cycle(
            cbcl_files[c], headers[c], filter_sums, n_tiles, n_reads
        )

        new_header = make_header(headers[c], new_tile_offsets)

        # make new file
        with (data_path / f"C{cycles_to_use[c]}.1" / "L001_1.cbcl").open("wb") as out:
            out.write(new_header)
            out.write(new_tiles)


if __name__ == "__main__":
    main()
