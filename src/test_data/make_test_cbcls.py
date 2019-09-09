#!/usr/bin/env python

import gzip
import pathlib
import struct

import xml.etree.ElementTree as et

import numpy as np

import click

from seqbot.bcl2fastr import NovaSeqRun, CBCLHeader, get_cycle


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
        key=lambda ffn: ffn.name.split("_")[-1],
    )

    for filter_file in filter_list[:n_tiles]:
        with filter_file.open("rb") as f:
            a, b, _ = struct.unpack("<III", f.read(12))
            n_filters = f.read(n_reads)

        with (data_path / filter_file.name).open("wb") as out:
            out.write(struct.pack("<III", a, b, n_reads))
            out.write(n_filters)


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
        header.bins.tobytes(),
        struct.pack("<I", new_tile_offsets.shape[0]),
        new_tile_offsets.tobytes(),
        struct.pack("B", header.non_PF_clusters_excluded),
    ]

    return b"".join(new_header)


def subset_cycle(
    cbcl_file: pathlib.Path, cbcl_header: CBCLHeader, n_tiles: int, n_reads: int
):
    """For a single CBCL file, read the first n_reads from the first n_tiles
    and return the gzip-compressed blocks and tile offsets
    """

    with cbcl_file.open("rb") as f:
        f.seek(cbcl_header.size)

        byte_arrays = []
        for i in range(n_tiles):
            byte_arrays.append(gzip.decompress(f.read(cbcl_header.tile_offsets[i, 3])))

    new_tiles = [gzip.compress(ba[: n_reads // 2]) for ba in byte_arrays]

    new_tile_offsets = np.array(
        [
            (cbcl_header.tile_offsets[i, 0], n_reads, n_reads // 2, len(new_tiles[i]))
            for i in range(n_tiles)
        ]
    )

    return b"".join(new_tiles), new_tile_offsets


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

    indexes = [(i1, i2) for i1, i2, index in run_info if index]
    ix_start = min(i1 for i1, i2 in indexes)
    ix_end = max(i2 for i1, i2 in indexes)

    # use n_bases / 2 from start so that we have some unfiltered files
    cycles_to_use = dict(zip(range(1, n_bases // 2 + 1), range(1, n_bases // 2 + 1)))

    # then add the indexes and another n_bases / 2 from second read
    for c in range(ix_start, ix_end + n_bases // 2 + 1):
        cycles_to_use[c] = len(cycles_to_use)

    cbcl_files = {
        get_cycle(cbcl_file): cbcl_file
        for cbcl_file in (
            run_path / "Data" / "Intensities" / "BaseCalls" / "L001"
        ).glob(f"C*.1/L001_1.cbcl")
        if get_cycle(cbcl_file) in cycles_to_use
    }

    _, _, headers, tiles = NovaSeqRun.read_headers(1, 1, cbcl_files)

    # write the new RunInfo file
    make_run_info(run_path, output_path, n_tiles, n_bases)

    data_path = output_path / "Data" / "Intensities" / "BaseCalls" / "L001"
    data_path.mkdir(parents=True, exist_ok=True)

    subset_locs(run_path, output_path, n_reads)

    subset_filters(run_path, output_path, n_tiles, n_reads)

    # write cycles
    for c in cycles_to_use:
        (data_path / f"C{cycles_to_use[c]}.1").mkdir(exist_ok=True)

        # take from lane 1 part 1
        new_tiles, new_tile_offsets = subset_cycle(
            cbcl_files[c], headers[c], n_tiles, n_reads
        )

        new_header = make_header(headers[c], new_tile_offsets)

        # make new file
        with (data_path / f"C{cycles_to_use[c]}.1" / "L001_1.cbcl").open("wb") as out:
            out.write(new_header)
            out.write(new_tiles)


if __name__ == "__main__":
    main()
