"""
    Reference Sequence Expansion Tool | RPINerd, 03/28/24

    With an alignment of sequences sometime the reference sequence is missing chuncks when compared
    to different variants. This tool will take in a fasta alignment file and treating the first entry as
    the "primary" it will create fastas for each segment (of configureable length) that is present
    in one of the other variant sequences.

    General Idea:

    RefSeq      ACCCACG-----TCCATCAA---GAAGTTCGCGTC------GATGATGCAGTGTGCTAGCTCGACACAGC
    VarSeq 1    ACCCACGATACGTCCATCAA---GAAGTTCGCGTC------GATGATGCAGTGT------TCGACACAGC
    VarSeq 2    ACCCACG-----TCCATCAATTCGAAGTTCGCGTC------GATGATGCAGTGT------TCGACACAGC
    VarSeq 3    ACCCACG-----TCCATCAATTCGAAGTTCGCGTCGGTGACGATGATGCAGTGT------TCGACACAGC
                        Exp1        Exp2            Exp3               N/A

    >Exp1
    ATACG
    >Exp2
    TTC
    >Exp3
    GGTGAC

"""

import argparse
import os
import sys


def parse_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta alignment file. First line must be reference!", default="test/ref_filler.fasta")
    args = parser.parse_args()
    return args


def validate_falign(file) -> str:
    """"""
    
    message = ""
    
    #TODO make sure there is more than 1 sequence
    #TODO input must be an alignment, i.e. have '-' characters

    return message


def parse_fasta(file: str) -> dict[str, set[str]]:
    """"""

    tracks = {}
    with open(file, "r") as alignment:
        lines = alignment.readlines()
        for id, seq in zip(lines[::2], lines[1::2]):
            
def main(args) -> None:

    # TODO validate fasta file instead of assuming input
    invalid = validate_falign(args.fasta)
    if invalid:
        raise f"Invalid fasta alignment input: {invalid}"
    
    # Extract the individual sequences as sets
    tracks = parse_fasta(args.fasta)

    return


if __name__ == "__main__":
    args = parse_args()
    main(args)
