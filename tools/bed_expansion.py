"""
Bed Expansion Tool | RPINerd, 01/21/26

Expand BED regions to a minimum length.

E.g. chr1:400-500 -> chr1:350-550
    100bp long region, the user requests a minimum length of 200bp:
    This program will expand the region to 150bp on either side (if possible) to reach the minimum length.
    If the region is already at least the minimum length, it will be left unchanged.
"""

import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        "%(asctime)s | %(levelname)-7s | %(module)-10s | %(lineno)-4d | %(message)s"
    ),
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def arg_parse() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Input BED file to expand"
    )
    parser.add_argument(
        "-n", "--num", type=int, help="Expand all regions to <number> minimum length"
    )
    args = parser.parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"{args.input.as_posix()} was not found!")

    return args


def expand(start: int, stop: int, target: int) -> tuple[int, int]:
    """Expand a region to a minimum length."""
    total_needed = target - abs(stop - start)
    end_split = total_needed // 2
    new_start = start - end_split
    new_stop = stop + end_split
    if abs(new_stop - new_start) < target:
        new_stop += 1
    return new_start, new_stop


def main(input: Path, minimum: int) -> None:
    """"""
    expanded = input.parent / f"{input.name.split('.')[0]}_expanded.bed"
    with (
        input.open("r", encoding="utf-8") as in_file,
        expanded.open("w", encoding="utf-8") as ex_file,
    ):
        for line in in_file:
            cols = line.strip().split()
            id = str(cols[0])
            start, stop = int(cols[1]), int(cols[2])
            if abs(stop - start) < minimum:
                ex_start, ex_stop = expand(start, stop, minimum)
                ex_file.write(f"{id}\t{ex_start}\t{ex_stop}\n")
            else:
                ex_file.write(line)


if __name__ == "__main__":
    args = arg_parse()
    main(args.input, args.num)
