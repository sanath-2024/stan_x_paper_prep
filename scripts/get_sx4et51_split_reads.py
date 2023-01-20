import sys
from dataclasses import dataclass
from typing import List

result_folder_name = sys.argv[1]


@dataclass
class Tgt:
    name: str
    te_name: str
    chrom: str
    pos: int


target = Tgt("SX4Et51", "copia", "2R", 9_237_984)

# number of bp on each side of the insertion that we allow the transposon to be
cutoff = 10_000

# list of relevant reads
reads: List[str] = []

with open(f"output/te_mapper/{result_folder_name}/genome_aligned.sam", "r") as in_file:
    # readlines: trailing newline is included
    for line in in_file.readlines():
        # skip header lines
        if line.startswith("@SQ") or line.startswith("@PG"):
            continue

        # extract useful info
        tabs = line.split("\t")
        te_name = tabs[0].split("|")[1]
        chrom = tabs[2]
        genome_pos = int(tabs[3])

        # only add insertions that correspond to the same element
        if (
            te_name.startswith(target.te_name)
            and chrom == target.chrom
            and abs(genome_pos - target.pos) <= cutoff
        ):
            reads.append(line)

print(
    f"found {len(reads)} reads for SX4Et51 in file output/te_mapper/{result_folder_name}/genome_aligned.sam"
)

with open(f"output/te_mapper/{result_folder_name}/sx4et51_reads.txt", "w") as out_file:
    for read in reads:
        out_file.write(read)