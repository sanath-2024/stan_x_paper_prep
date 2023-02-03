import re
import sys
from dataclasses import dataclass
from typing import List, Tuple

from PIL import Image, ImageDraw

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

IMAGE_HEIGHT_MARGIN_EACH_SIDE = 100
IMAGE_WIDTH_MARGIN_EACH_SIDE = 200
NT_PER_PIXEL = 1
BAR_HEIGHT = 5
BAR_SPACING = 3
INSERTION_COLOR = (0, 0, 255)
TE_PART_COLOR = (0, 255, 0)
GENOME_PART_COLOR = (255, 0, 0)

SX4ET51_TRANSPOSON_NAME = "copia#LTR/Copia"
SX4ET51_TRANSPOSON_LENGTH = 5143


def draw_single_insertion(
    alignments: List[UpstreamTEAlignment],
    draw: ImageDraw,
    lower_left_pos: Tuple[int, int],
):
    insertion_size = SX4ET51_TRANSPOSON_LENGTH

    # sort the upstream and downstream reads to get the V shape when we iterate through them
    alignments.sort(key=lambda x: x.te_range[0])

    # draw transposon
    draw.rectangle(
        [
            (lower_left_pos[0], lower_left_pos[1]),
            (lower_left_pos[0] + insertion_size, lower_left_pos[1] - BAR_HEIGHT),
        ],
        fill=INSERTION_COLOR,
    )

    # draw upstream reads
    cursor = (lower_left_pos[0], lower_left_pos[1])
    for read in alignments:
        cursor = (cursor[0], cursor[1] - BAR_HEIGHT - BAR_SPACING)
        draw.rectangle(
            [
                (
                    cursor[0] - (read.genome_range[1] - read.genome_range[0]) + 1,
                    cursor[1],
                ),
                (cursor[0], cursor[1] - BAR_HEIGHT),
            ],
            fill=GENOME_PART_COLOR,
        )
        draw.rectangle(
            [
                (cursor[0], cursor[1]),
                (
                    cursor[0] + read.te_range[1] - read.te_range[0] + 1,
                    cursor[1] - BAR_HEIGHT,
                ),
            ],
            fill=TE_PART_COLOR,
        )


def main():
    insertion_size = SX4ET51_TRANSPOSON_LENGTH
    image_width = insertion_size + 2 * IMAGE_WIDTH_MARGIN_EACH_SIDE
    image_height = (
        BAR_HEIGHT
        + len(alignments) * (BAR_HEIGHT + BAR_SPACING)
        + 2 * IMAGE_HEIGHT_MARGIN_EACH_SIDE
    )
    image = Image.new("RGB", (image_width, image_height), (255, 255, 255))
    draw_single_insertion(
        alignments,
        ImageDraw.Draw(image),
        (
            IMAGE_WIDTH_MARGIN_EACH_SIDE,
            image_height - IMAGE_HEIGHT_MARGIN_EACH_SIDE,
        ),
    )
    image.save(f"output/te_mapper/{result_folder_name}/figure_sx4et51_upstream.png")


if __name__ == "__main__":
    main()
