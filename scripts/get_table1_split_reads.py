import json
import sys
from dataclasses import dataclass
from itertools import zip_longest
from typing import List

result_folder_name = sys.argv[1]


@dataclass
class Tgt:
    name: str
    te_name: str
    chrom: str
    pos: int


# list of SX-4 insertions within a natural TE,
# found through inverse PCR (taken from Google Sheet)
targets: List[Tgt] = [
    Tgt("SX4Ch7", "1360", "2L", 12_004_570),
    Tgt("SX4Aq839", "1360", "2L", 16_727_570),
    Tgt("SX4Lv807", "invader1", "2R", 6_622_465),
    Tgt("SX4Et51", "copia", "2R", 9_237_984),
    Tgt("SX4Et8", "HMS-Beagle", "2R", 15_951_007),
    Tgt("SX4Et49", "opus", "3L", 17_918_916),
    Tgt("SX4Lv831", "Juan", "3R", 4_411_749),
    Tgt("SX4Lv816", "1360", "3R", 4_610_657),
    Tgt("SX4Lv811", "F-element", "3R", 4_753_706),
    Tgt("SX4Co882", "mdg3", "3R", 5_073_316),
    Tgt("SX4ECPS11", "invader4", "3R", 16_189_617),
]

# number of bp on each side of the insertion that we allow the transposon to be
cutoff = 1_000


@dataclass
class SplitReadRanges:
    te_range: List[int]
    genome_range: List[int]


@dataclass
class Result:
    ref: bool
    name: str
    chrom: str
    upstream_pos: int
    downstream_pos: int
    orientation: str
    upstream_reads: List[SplitReadRanges]
    downstream_reads: List[SplitReadRanges]


with open(
    f"output/te_mapper/{result_folder_name}/te_mapper_output.json", "r"
) as in_file:
    raw = json.load(in_file)
    reads = []
    for chrom in raw:
        for read in chrom["non_reference"]:
            reads.append(Result(False, **read))
        for read in chrom["reference"]:
            reads.append(Result(True, **read))
    for read in reads:
        read.upstream_reads = [SplitReadRanges(**x) for x in read.upstream_reads]
        read.downstream_reads = [SplitReadRanges(**x) for x in read.downstream_reads]


def candidate(tgt: Tgt, res: Result) -> bool:
    """Is the result a candidate to be the target?"""

    # check that it is the same transposon
    if not res.name.startswith(tgt.te_name):
        return False

    # check that it is on the same chromosome
    if res.chrom != tgt.chrom:
        return False

    # check the position
    if res.ref:
        return tgt.pos in range(res.upstream_pos, res.downstream_pos + 1)
    else:
        return (
            abs(read.upstream_pos - tgt.pos) <= cutoff
            or abs(read.downstream_pos - tgt.pos) <= cutoff
        )


filtered: List[Result] = []
for tgt in targets:
    # for some reason, list comprehension doesn't work here
    tentative = []
    for read in reads:
        if candidate(tgt, read):
            tentative.append(read)
    if len(tentative) > 1:
        raise RuntimeError(
            f"tentative insertion list has more than 1 element (target: {tgt})"
        )
    elif len(tentative) == 1:
        filtered.append((tgt, tentative[0]))

print([f[0].name for f in filtered])

with open(
    f"output/te_mapper/{result_folder_name}/split_reads_for_tgt.csv", "w"
) as out_file:
    out_file.write(
        "Name,"
        "TE Name,"
        "Chromosome,"
        "Upstream Position,"
        "Downstream Position,"
        "Orientation,"
        "Reference?,"
        "Upstream Reads (TE Range),,"
        "Upstream Reads (Genome Range),,"
        "Downstream Reads (TE Range),,"
        "Downstream Reads (Genome Range),\n"
    )
    for tgt, result in filtered:
        out_file.write(
            f"{tgt.name},"
            f"{result.name},"
            f"{result.chrom},"
            f"{result.upstream_pos},"
            f"{result.downstream_pos},"
            f"{'+/+' if result.orientation == 'PlusPlus' else '+/-'},"
            f"{'reference' if result.ref else 'non-reference'},"
            ",,,,,,,\n"
        )
        reads = zip_longest(result.upstream_reads, result.downstream_reads)
        for upstream, downstream in reads:
            out_file.write(",,,,,,,")
            if upstream is None:
                out_file.write(",,,,")
            else:
                out_file.write(
                    f"{upstream.te_range[0]},{upstream.te_range[1]},{upstream.genome_range[0]},{upstream.genome_range[1]},"
                )
            if downstream is None:
                out_file.write(",,,\n")
            else:
                out_file.write(
                    f"{downstream.te_range[0]},{downstream.te_range[1]},{downstream.genome_range[0]},{downstream.genome_range[1]}\n"
                )
