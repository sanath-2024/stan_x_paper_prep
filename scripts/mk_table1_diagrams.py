import json
from dataclasses import dataclass
from typing import List, Optional, Tuple

from PIL import Image, ImageDraw


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


@dataclass
class Target:
    name: str
    upstream_pos: int
    length: Optional[int] = None
    result: Optional[Result] = None


targets = [
    Target("SX4Ch7", 12_004_481, 3409),
    Target("SX4Et8", 15_951_006, 7062),
    Target("SX4Et49", 17_918_427),
    Target("SX4Lv816", 4_609_605),
]


def get_results() -> List[Result]:
    with open("output/te_mapper/1_R1_001/te_mapper_output.json", "r") as in_file:
        raw = json.load(in_file)
        reads = []
        for chrom in raw:
            for read in chrom["non_reference"]:
                reads.append(Result(False, **read))
            for read in chrom["reference"]:
                reads.append(Result(True, **read))
        for read in reads:
            read.upstream_reads = [SplitReadRanges(**x) for x in read.upstream_reads]
            read.downstream_reads = [
                SplitReadRanges(**x) for x in read.downstream_reads
            ]
        return reads


def filter_results(results: List[Result]):
    for target in targets:
        found = False
        for result in results:
            if result.upstream_pos == target.upstream_pos:
                found = True
                target.result = result
                break
        if not found:
            raise RuntimeError(f"insertion {target.name} not found")


IMAGE_HEIGHT_MARGIN_EACH_SIDE = 100
IMAGE_WIDTH_MARGIN_EACH_SIDE = 200
NT_PER_PIXEL = 1
BAR_HEIGHT = 5
BAR_SPACING = 3
INSERTION_COLOR = (0, 0, 255)
TE_PART_COLOR = (0, 255, 0)
GENOME_PART_COLOR = (255, 0, 0)


def draw_single_insertion(
    target: Target, draw: ImageDraw, lower_left_pos: Tuple[int, int]
):
    insertion_size = (
        target.length
        if target.length is not None
        else target.result.downstream_pos - target.result.upstream_pos + 1
    )

    insertion = target.result

    # sort the upstream and downstream reads to get the V shape when we iterate through them
    insertion.upstream_reads.sort(key=lambda x: x.te_range[0])
    insertion.downstream_reads.sort(key=lambda y: y.te_range[1], reverse=True)

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
    for read in insertion.upstream_reads:
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

    # draw downstream reads
    cursor = (lower_left_pos[0] + insertion_size, lower_left_pos[1])
    for read in insertion.downstream_reads:
        cursor = (cursor[0], cursor[1] - BAR_HEIGHT - BAR_SPACING)
        draw.rectangle(
            [
                (
                    cursor[0] - (read.te_range[1] - read.te_range[0]) + 1,
                    cursor[1],
                ),
                (cursor[0], cursor[1] - BAR_HEIGHT),
            ],
            fill=TE_PART_COLOR,
        )
        draw.rectangle(
            [
                (cursor[0], cursor[1]),
                (
                    cursor[0] + read.genome_range[1] - read.genome_range[0] + 1,
                    cursor[1] - BAR_HEIGHT,
                ),
            ],
            fill=GENOME_PART_COLOR,
        )


def main():
    filter_results(get_results())
    for target in targets:
        result = target.result
        insertion_size = (
            target.length
            if target.length is not None
            else target.result.downstream_pos - target.result.upstream_pos + 1
        )
        image_width = insertion_size + 2 * IMAGE_WIDTH_MARGIN_EACH_SIDE
        image_height = (
            BAR_HEIGHT
            + max(len(result.upstream_reads), len(result.downstream_reads))
            * (BAR_HEIGHT + BAR_SPACING)
            + 2 * IMAGE_HEIGHT_MARGIN_EACH_SIDE
        )
        image = Image.new("RGB", (image_width, image_height), (255, 255, 255))
        draw_single_insertion(
            target,
            ImageDraw.Draw(image),
            (
                IMAGE_WIDTH_MARGIN_EACH_SIDE,
                image_height - IMAGE_HEIGHT_MARGIN_EACH_SIDE,
            ),
        )
        image.save(f"output/te_mapper/1_R1_001/figure_{target.name}.png")


if __name__ == "__main__":
    main()
