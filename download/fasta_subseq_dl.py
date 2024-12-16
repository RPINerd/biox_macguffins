"""
    Fasta Subsequence Download Tool | RPINerd, 11/26/24

    Given a file with Accession IDs and start/stop coordinates, download the fasta sequence for each entry.
    Entry format is as follows:
        CP060352.1:280230-280675
"""

import argparse
import sys
from pathlib import Path
from time import sleep

from Bio import Entrez


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Input file with Accession IDs and start/stop coordinates")
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    """Main function"""
    if not Path.exists(args.input):
        print(f"Error: {args.input} does not exist.")
        sys.exit(1)

    # TODO refactor to submit all IDs as single request
    Entrez.email = "rpinerd@gmail.com"
    with Path.open(args.input) as id_list:
        sleep(1)
        for line in id_list:
            enst, coords = line.strip().split(":")
            start, stop = coords.split("-")
            handle = Entrez.efetch(db="nucleotide", id=enst, rettype="fasta")
            record = handle.read()
            handle.close()
            lines = record.split("\n")
            header = lines[0]
            sequence = "".join(lines[1:])

            with Path.open(f"{enst.split(".")[0]}.fasta", "w") as out:
                out.write(f"{header}\n{sequence[int(start):int(stop)]}\n")


if __name__ == "__main__":
    args = parse_args()
    main(args)
