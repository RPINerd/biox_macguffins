"""
    Degenerate Bases | RPINerd, 03/06/24

    Given a directory, this script will search all fasta files for degenerate bases (non A, T, C or G)
    and report the number found.
"""

import os
import sys

from Bio import SeqIO


def main(directory) -> None:

    # Step through the directory and search for fasta files
    for file in os.listdir(directory):
        if file.endswith(".fna") or file.endswith(".fa") or file.endswith(".fasta"):
            # Parse the fasta file using biopython
            records = list(SeqIO.parse(os.path.join(directory, file), "fasta"))
            for record in records:
                # Count any base that is not A T C or G
                ambig = 0
                for base in record.seq:
                    if base not in ["A", "T", "C", "G"]:
                        ambig += 1

                print(f"{file:25}{record.id:20}{len(record.seq):10}{ambig:8}")
    return


if __name__ == "__main__":
    main(sys.argv[1])
