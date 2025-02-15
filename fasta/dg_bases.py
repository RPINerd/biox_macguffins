"""
    Degenerate Bases | RPINerd, 03/06/24

    Given a directory, this script will search all fasta files for degenerate bases (non A, T, C or G)
    and report the number found.
"""

import os
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord


def main(directory: str) -> None:
    """
        Main function

    :param directory: Directory containing fasta files

    :return: None
    """
    # Step through the directory and search for fasta files
    for file in os.listdir(directory):
        if file.endswith(".fna") or file.endswith(".fa") or file.endswith(".fasta"):
            # Parse the fasta file using biopython
            records: list[SeqRecord] = list(SeqIO.parse(Path(directory) / file / "fasta"))
            for record in records:
                # Count any base that is not A T C or G
                ambig = 0
                for base in record.seq:
                    if base not in {"A", "T", "C", "G"}:
                        ambig += 1

                print(f"{file:25}{record.id:20}{len(record.seq):10}{ambig:8}")


if __name__ == "__main__":
    main(sys.argv[1])
