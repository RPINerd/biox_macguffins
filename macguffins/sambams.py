"""Functions related to manipulation of sam/bam files"""

import re
import sys
from pathlib import Path


def cigar_filter(infile: Path, filter_value: int = 50, outfile: Path | None = None) -> None:
    """
    Barebones function to filter out reads based on a single CIGAR characteristic

    Args:
        infile (Path): Path to the input SAM file
        filter_value (int): Minimum value for the CIGAR filter
        outfile (Path): Path to the output file. If None, output is written to stdout
    """
    if outfile:
        out = Path.open(outfile, "w")
    else:
        out = sys.stdout

    filtered_reads: list[list[str]] = []
    with Path.open(infile, "r") as f:
        for line in f:

            # Skip the header info lines
            if line.startswith("@"):
                continue

            cols = line.split("\t")
            filter = re.match(r"([0-9]+)M", cols[5])  # column 5 is the cigar
            if filter and int(filter[1]) >= filter_value:
                filtered_reads.append(cols)
                out.write("\t".join([cols[2], cols[3], cols[5]]) + "\n")
