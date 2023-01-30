import re
import sys
from dataclasses import dataclass
from typing import List

result_folder_name = sys.argv[1]


@dataclass
class Tgt:
    name: str
    chrom: str
    pos: int


target = Tgt("SX4Et51", "2R", 9_237_984)

# number of bp on each side of the insertion that we allow the transposon to be
cutoff = 500

# list of relevant reads
reads: List[str] = []

with open(f"output/bwa_genome/{result_folder_name}/genome_aligned.sam", "r") as in_file:
    # readlines: trailing newline is included
    for line in in_file.readlines():
        # skip header lines
        if line.startswith("@SQ") or line.startswith("@PG"):
            continue

        # extract useful info
        tabs = line.split("\t")
        chrom = tabs[2]
        genome_pos = int(tabs[3])
        cigar_string = tabs[5]
        expr = re.compile(r"\d+(S|H)\d+M")

        # only get reads that aren't 150M
        if (
            chrom == target.chrom
            and abs(genome_pos - target.pos) <= cutoff
            and expr.fullmatch(cigar_string)
        ):
            reads.append(line)

print(
    f"found {len(reads)} reads for SX4Et51 in file output/bwa_genome/{result_folder_name}/genome_aligned.sam"
)

with open(f"output/bwa_genome/{result_folder_name}/sx4et51_reads.txt", "w") as out_file:
    for read in reads:
        out_file.write(read)
