import re
import sys
from dataclasses import dataclass
from typing import List, Tuple

result_folder_name = sys.argv[1]


@dataclass
class UpstreamTEAlignment:
    te_range: Tuple[int, int]
    genome_range: Tuple[int, int]


with open(f"output/te_mapper/{result_folder_name}/sx4et51_reads.txt", "r") as in_file:
    alignments: List[UpstreamTEAlignment] = []
    for line in in_file.readlines():
        fields = line.split("\t")
        cigar_string = fields[5]
        expr = re.compile(r"(\d+)M(\d+)(S|H)")
        matches = expr.fullmatch(cigar_string)
        if matches is None:
            continue
        genome_range_size, te_range_size = matches.group(1, 2)
        te_range_size = int(te_range_size)
        genome_range_size = int(genome_range_size)
        alignments.append(
            UpstreamTEAlignment(
                (genome_range_size + 1, genome_range_size + te_range_size),
                (1, genome_range_size),
            )
        )

with open(f"output/te_mapper/{result_folder_name}/sx4et51_reads.csv", "w") as out_file:
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
    out_file.write(
        "SX4Et51,"
        "copia#LTR/Copia,"
        "2R,"
        "9237984,"
        ","
        "+/+,"
        "non-reference,"
        ",,,,,,,\n"
    )
    for alignment in alignments:
        out_file.write(",,,,,,,")
        out_file.write(
            f"{alignment.te_range[0]},{alignment.te_range[1]},{alignment.genome_range[0]},{alignment.genome_range[1]},"
        )
        out_file.write(",,,\n")
