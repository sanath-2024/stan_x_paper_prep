import json
import re
from dataclasses import dataclass
from typing import Dict, List, Optional

result_folders = {
    "output/te_mapper/1_R1_001": "1_R1_001.fastq",
    "output/te_mapper/1_R2_001": "1_R2_001.fastq",
    "output/te_mapper/2_R1_001": "2_R1_001.fastq",
    "output/te_mapper/2_R2_001": "2_R2_001.fastq",
    "output/te_mapper/7_R1_001": "7_R1_001.fastq",
    "output/te_mapper/7_R2_001": "7_R2_001.fastq",
    "output/te_mapper/8_R1_001": "8_R1_001.fastq",
    "output/te_mapper/8_R2_001": "8_R2_001.fastq",
}

result_folders_inv = {
    "1_R1_001.fastq": "output/te_mapper/1_R1_001",
    "1_R2_001.fastq": "output/te_mapper/1_R2_001",
    "2_R1_001.fastq": "output/te_mapper/2_R1_001",
    "2_R2_001.fastq": "output/te_mapper/2_R2_001",
    "7_R1_001.fastq": "output/te_mapper/7_R1_001",
    "7_R2_001.fastq": "output/te_mapper/7_R2_001",
    "8_R1_001.fastq": "output/te_mapper/8_R1_001",
    "8_R2_001.fastq": "output/te_mapper/8_R2_001",
}


@dataclass
class Tgt:
    name: str
    te_name: str
    chrom: str
    pos: int


# list of relevant SX-4 insertions,
# confirmed through inverse PCR
targets: List[Tgt] = [
    Tgt("SX4Ch7", "1360#DNA/P", "2L", 12_004_570),
    Tgt("SX4Aq839", "1360#DNA/P", "2L", 16_727_570),
    Tgt("SX4Et51", "copia#LTR/Copia", "2R", 9_237_984),
    Tgt("SX4Et8", "HMS-Beagle#LTR/Gypsy", "2R", 15_951_007),
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
    target_name: Optional[str] = None
    te_name: Optional[str] = None


results: Dict[str, List[Result]] = {}

for result_folder in result_folders:
    results[result_folders[result_folder]] = []
    with open(f"{result_folder}/te_mapper_output.json", "r") as in_file:
        raw = json.load(in_file)
        for chrom in raw:
            for result in chrom["non_reference"]:
                results[result_folders[result_folder]].append(Result(False, **result))
            for result in chrom["reference"]:
                results[result_folders[result_folder]].append(Result(True, **result))
        for result in results[result_folders[result_folder]]:
            result.upstream_reads = [
                SplitReadRanges(**x) for x in result.upstream_reads
            ]
            result.downstream_reads = [
                SplitReadRanges(**x) for x in result.downstream_reads
            ]


def candidate(tgt: Tgt, res: Result) -> bool:
    """Is the result a candidate to be the target?"""

    # check that it is the same transposon
    if res.name != tgt.te_name:
        return False

    # check that it is on the same chromosome
    if res.chrom != tgt.chrom:
        return False

    # check the position
    if res.ref:
        return tgt.pos in range(res.upstream_pos, res.downstream_pos + 1)
    else:
        return (
            abs(res.upstream_pos - tgt.pos) <= cutoff
            or abs(res.downstream_pos - tgt.pos) <= cutoff
        )


filtered: Dict[str, List[Result]] = {}

for filename in results:
    filtered[filename] = []
    for tgt in targets:
        tentative = [result for result in results[filename] if candidate(tgt, result)]
        if len(tentative) > 1:
            raise RuntimeError(
                f"tentative insertion list has more than 1 element (target: {tgt})"
            )
        elif len(tentative) == 1:
            tentative[0].target_name = tgt.name
            tentative[0].te_name = tgt.te_name
            filtered[filename].append(tentative[0])


@dataclass
class Row:
    fastq_filename: str
    insertion_line: str
    natural_te: str
    sequence: str
    upstream: bool


rows: List[Row] = []

upstream_regex = re.compile(r"(\d+)M(\d+)(S|H)")
downstream_regex = re.compile(f"(\d+)(S|H)(\d+)M")

for filename in filtered:
    print(f"processing genome alignments for file {filename} ...")
    with open(
        f"{result_folders_inv[filename]}/genome_aligned.sam", "r"
    ) as genome_aligned_file:
        alignments = genome_aligned_file.read().splitlines()
        alignments = [line for line in alignments if not line.startswith("@")]
    for line in alignments:
        sam_fields = line.split("\t")

        # parse out genome alignment information
        cigar_string = sam_fields[5]
        upstream_match = upstream_regex.fullmatch(cigar_string)
        downstream_match = downstream_regex.fullmatch(cigar_string)
        if upstream_match is None and downstream_match is None:
            continue

        chrom = sam_fields[2]
        pos = int(sam_fields[3])
        seq = sam_fields[9]
        if len(seq) != 150:
            continue

        if upstream_match is not None:
            upstream = True
            genome_match_size, genome_clip_size = upstream_match.group(1, 2)
        else:
            upstream = False
            genome_match_size, genome_clip_size = downstream_match.group(3, 1)
        genome_match_size = int(genome_match_size)
        genome_clip_size = int(genome_clip_size)

        # parse out TE alignment information
        te_alignment_info = sam_fields[0].split("|")
        te_name = te_alignment_info[1]
        te_match_size = int(te_alignment_info[2])
        te_clip_size = int(te_alignment_info[3])

        if genome_match_size != te_clip_size or genome_clip_size != te_match_size:
            continue

        for result in filtered[filename]:
            if te_name != result.te_name:
                continue
            if chrom != result.chrom:
                continue
            if (
                abs(result.upstream_pos - pos) > cutoff
                and abs(result.downstream_pos - pos) > cutoff
            ):
                continue
            rows.append(Row(filename, result.target_name, te_name, seq, upstream))

with open("output/te_mapper/split_read_sequence_table.csv", "w") as out_file:
    out_file.write(
        "FASTQ File Name,Insertion Line Name,Natural TE,Match End,Sequence\n"
    )
    for row in rows:
        out_file.write(
            f"{row.fastq_filename},{row.insertion_line},{row.natural_te},{'Upstream' if row.upstream else 'Downstream'},{row.sequence}\n"
        )
