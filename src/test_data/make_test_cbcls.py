#!/usr/bin/env python

import gzip
import pathlib
import struct

import xml.etree.ElementTree as et

import numpy as np

import click

from seqbot.bcl2fastr import NovaSeqRun, CBCLHeader


def make_run_info_file(
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


def make_header(header: CBCLHeader, new_tile_offsets: np.ndarray):
    header_size = 17 + 4 * header.bins.size + 4 * new_tile_offsets.size

    new_header = [
        struct.pack(
            "<HIBBI", header.version, header_size, header.bits_per_basecall,
            header.bits_per_qscore, header.number_of_bins
        ),
        header.bins.tobytes(),
        struct.pack("<I", new_tile_offsets.shape[0]),
        new_tile_offsets.tobytes(),
        struct.pack("B", header.non_PF_clusters_excluded)
    ]

    return b"".join(new_header)


def subset_cycle(
    novaseq_run: NovaSeqRun,
    cycle: int,
    lane: int,
    part: int,
    n_tiles: int,
    n_reads: int,
):
    """For a single CBCL file, read the first n_reads from the first n_tiles
    and return the gzip-compressed blocks

    :param novaseq_run:
    :param cycle:
    :param lane:
    :param part:
    :param n_tiles:
    :param n_reads:
    :return:
    """

    c_h = novaseq_run.headers[lane, part][cycle]

    with novaseq_run.cbcl_files[lane, part][cycle].open("rb") as f:
        f.seek(c_h.size)

        byte_arrays = []
        for i in range(n_tiles):
            byte_arrays.append(gzip.decompress(f.read(c_h.tile_offsets[i, 3])))

    new_tiles = [gzip.compress(ba[:n_reads // 2]) for ba in byte_arrays]

    new_tile_offsets = np.array(
        [
            (c_h.tile_offsets[i, 0], n_reads, n_reads // 2, len(new_tiles[i]))
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
    "--output_dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("--n_tiles", required=True, type=int)
@click.option("--n_reads", required=True, type=int)
@click.option("--n_bases", required=True, type=int)
def main(
    run_path: str,
    output_dir: str,
    n_tiles: int,
    n_reads: int,
    n_bases: int,
):
    assert n_reads % 2 == 0

    run_path = pathlib.Path(run_path)
    output_dir = pathlib.Path(output_dir)

    novaseq_run = NovaSeqRun(run_path)

    indexes = [(i1,i2) for i1, i2, index in novaseq_run.run_info if index]
    ix_start = min(i1 for i1, i2 in indexes)
    ix_end = max(i2 for i1, i2 in indexes)

    # use n_bases / 2 from start so that we have some unfiltered files
    cycles_to_use = dict(zip(range(1, n_bases // 2 + 1), range(1, n_bases // 2 + 1)))

    # then add the indexes and another n_bases / 2 from second read
    for c in range(ix_start, ix_end + n_bases // 2 + 1):
        cycles_to_use[c] = len(cycles_to_use)

    data_path = output_dir / "Data" / "Intensities" / "BaseCalls" / "L001"
    data_path.mkdir(parents=True, exist_ok=True)

    # write the new RunInfo file
    make_run_info_file(run_path, output_dir, n_tiles, n_bases)

    # read raw locs
    with (run_path / "Data" / "Intensities" / "s.locs").open("rb") as f:
        a, b, _ = struct.unpack("<III", f.read(12))
        cbcl_locs = f.read(8 * n_reads)

    # write locs to new file
    with (output_dir / "Data" / "Intensities" / "s.locs").open("wb") as out:
        out.write(struct.pack("<III", a, b, n_reads))
        out.write(cbcl_locs)

    # write filter file
    for tile_i in novaseq_run.tiles[1, 1][:n_tiles]:
        with (data_path / f"s_1_{tile_i}.filter").open("wb") as out:
            out.write(struct.pack("<III", 0, 1, n_reads))
            out.write(novaseq_run.filters[1][tile_i][:n_reads].tobytes())

    # write cycles
    for c in cycles_to_use:
        (data_path / f"C{cycles_to_use[c]}.1").mkdir(exist_ok=True)
        # take from lane 1 part 1
        new_tiles, new_tile_offsets = subset_cycle(
            novaseq_run, c, 1, 1, n_tiles, n_reads
        )

        new_header = make_header(novaseq_run.headers[1, 1][c], new_tile_offsets)

        # make new file
        with (data_path / f"C{cycles_to_use[c]}.1" / "L001_1.cbcl").open("wb") as out:
            out.write(new_header)
            out.write(new_tiles)


if __name__ == "__main__":
    main()
