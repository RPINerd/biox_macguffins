"""
    BLAST Tools | RPINerd, 05/20/21

    Tools for manipulating/working with blast data
"""

import argparse
import re
from pathlib import Path


def self_ref_remove(in_lines: list[str]) -> list[str]:
    """
    Remove self hits from a blastn report file

    Args:
        in_lines (list[str]): List of lines from the blastn report file

    Returns:
        list[str]: List of lines with self hits removed
    """
    r_str = r"^(chr.+:[0-9]+\-[0-9]+)\t(chr.+:[0-9]+\-[0-9]+)"
    out = []

    total = len(in_lines)
    for prog, line in enumerate(in_lines):
        # Feedback
        current_prog = int((prog / total) * 100)
        if current_prog % 5 == 0:
            print("Removing self hits.. " + str(current_prog) + "%", end="\r")

        hit = re.match(r_str, line)
        if hit.group(1) == hit.group(2):
            next
        else:
            out.append(str(line))

    print("Removing self hits.. Complete!")
    return out


def pick_best(in_lines: list[str]) -> list[str]:
    out = []

    return out


def write_output(out_lines: list[str], out_file: str) -> None:
    """
    Given the output lines and file, write the output

    Args:
        out_lines (list[str]): List of lines to write to the output file
        out_file (str): Path to the output file

    Returns:
        None
    """
    with Path.open(out_file, "w") as out:
        for line in out_lines:
            out.write(str(line))


def main(args: argparse.Namespace) -> None:
    """"""
    data = []
    with Path.open(args.input) as file:
        for line in file:
            data.append(line)

    if args.trimselfhits:
        data = self_ref_remove(data)
    if args.pickbesthit:
        data = pick_best(data)

    write_output(data, args.output)


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    parser.add_argument("-o", "--output", help="Designates a user-defined output file", default="output.txt")
    parser.add_argument("--trimselfhits", action="store_true")
    parser.add_argument("--pickbesthit", action="store_true")
    args = parser.parse_args()

    main(args)
