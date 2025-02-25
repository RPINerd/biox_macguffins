"""
    CIGAR Filter | RPINerd, 12/09/24

    For now just a
"""

import re
import sys
from pathlib import Path


def main(filter_value: int = 50) -> None:
    """Main function"""
    out = Path.open(sys.argv[2], "w")

    for line in Path.open(sys.argv[1], "r"):
        if line.startswith("@"):
            continue
        cols = line.split("\t")
        # cols[5] is the cigar
        filter = re.match(r"([0-9]+)M", cols[5])
        if filter and int(filter[1]) >= filter_value:
            print(cols[2], cols[3], cols[5])
            out.write("\t".join([cols[2], cols[3], cols[5]]) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
