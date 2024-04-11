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

from simples import look_forward

# Minimum gap size in refseq to consider filling
REGION_SIZE = 150
# Bases up/down stream of patch to include
UPSTREAM_ANCHOR = 100
DOWNSTREAM_ANCHOR = 100


def parse_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        help="Input fasta alignment file. First line must be reference!",
        default="test/ref_filler.fasta",
    )
    parser.add_argument(
        "-g",
        "--gap",
        type=int,
        help="Minimum gap size to consider filling (Default = 150)",
        default=150,
        required=False,
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        help="Number of bases upstream to include in the patch fasta (Default = 100)",
        default=100,
        required=False,
    )
    parser.add_argument(
        "-d",
        "--downstream",
        type=int,
        help="Number of bases downstream to include in the patch fasta (Default = 100)",
        default=100,
        required=False,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output directory for expansion fasta files (Defaults to current directory)",
        default="./",
        required=False
    )
    args = parser.parse_args()
    return args


# TODO make sure there is more than 1 sequence
# TODO input must be an alignment, i.e. have '-' characters
def validate_falign(file) -> str:
    """
    Validate the input fasta alignment file

    :param file: Path to the fasta alignment file
    :rtype: str - Error message if invalid, else None
    """

    message = ""

    return message


def parse_fasta(file: str) -> dict[str, tuple[str]]:
    """
    Parse a fasta alignment file into a dictionary of sequences
    Sequences are stored as tuples of the individual bases

    :param str: file: Path to the fasta alignment file
    :rtype: dict[str, tuple[str]] - Dictionary of sequences with the ID as the key
    """

    tracks = {}
    with open(file, "r") as alignment:
        lines = alignment.readlines()

        ref = True
        for id, seq in zip(lines[::2], lines[1::2]):
            if ref:
                tracks["refseq"] = tuple(seq)
                ref = False
            else:
                tracks[id] = tuple(seq)

    return tracks


def parse_refseq(ref_seq: tuple[str]) -> list[tuple[int, int]]:
    """
    Parse the reference sequence to identify missing regions

    :param tuple[str]: ref_seq: Tuple of the reference sequence
    :rtype: list[tuple[int, int]] - List of tuples containing the start and end of missing regions
    """

    # TODO what if we walk through in steps of REGION_SIZE, since we know that anything contained within that scope would be too small to qualify anyways, then look_backward to find the start and look_forward for the end when landing on a gap

    # Unpack refseq and set up some variables
    sequence = ref_seq
    gaps = []
    gap_start = None
    gap_end = None
    nt_total = len(sequence)

    # -print(id, sequence)
    idx = 0
    while idx < nt_total:

        # -print(idx, sequence[idx])

        if sequence[idx] != "-":  # Nothing to see here
            idx += 1

        # A gap base is found at the current index
        else:

            # -print("gap base")
            # Verify gap length is valid
            minpos = idx + REGION_SIZE - 1
            # -print(f"minpos: {minpos}")
            # -print(f"subseq: {sequence[idx:minpos]}")
            # TODO could just infer size validity from the look_forward number, reduce duplicate lookup
            if not all(nt == "-" for nt in sequence[idx:minpos]):
                continue

            gap_start = idx
            # -print(f"start: {idx}, {sequence[idx]}")

            # Look forward until the next non-gap base is found
            gap_end = look_forward(sequence, minpos, "-")

            # Add this gap to the list
            gaps.append((gap_start, gap_end))

            # Advance the index to the next base of the seq following the gap
            idx = gap_end + 1

    return gaps


def main(args) -> None:

    fasta_file = args.fasta
    output_directory = args.output

    # TODO validate fasta file instead of assuming input
    invalid = validate_falign(fasta_file)
    if invalid:
        raise f"Invalid fasta alignment input: {invalid}"

    # Extract the individual sequences as sets
    tracks = parse_fasta(fasta_file)

    # Process the reference sequence to establish relevant gaps
    ref_seq = tracks.pop("refseq")
    missing_regions = parse_refseq(ref_seq)

    #! Debug
    print(f" Missing regions: {missing_regions}")

    # Process the other sequences to extract the missing regions
    algo = 2
    # - Option 1
    """
        for each missing region,
            get the midpoint index and
            for each other sequences in tracks
                check to see if it is a valid base or a gap (i.e. "-")

    """
    if algo == 1:
        for region in missing_regions:
            start, end = region
            midpoint = (start + end) // 2

            #! Debug
            print(midpoint)

            for id, seq in tracks.items():
                if seq[midpoint] != "-":
                    #! Debug
                    print(id, seq[midpoint])

    # - Option 2
    """
        for each missing region
            check that all positions of tracks[id] are not "-" within the region
    """
    if algo == 2:
        fasta_number = 1
        for region in missing_regions:
            start, end = region
            start_upstream = start - UPSTREAM_ANCHOR
            end_downstream = end + DOWNSTREAM_ANCHOR

            for id, seq in tracks.items():
                if all(nt != "-" for nt in seq[start:end]):
                    # TODO figure out how to best name these subsequences

                    #! Debug
                    print(f"Filler from: {id.strip()} with seq {"".join(seq[start:end])}")

                    # TODO try except for when the downstream/upstream are outside of the range, default to just the start or end
                    fasta_header = f">ExpansionSEQ {fasta_number}"
                    fasta_seq = "".join(seq[start_upstream:end_downstream])

                    #! Debug
                    print(f"Seq including anchors:\n{fasta_seq}")

                    with open(f"{output_directory}expansion{fasta_number}.fasta", "x") as expansion_file:
                        expansion_file.write(fasta_header + "\n")
                        expansion_file.write(fasta_seq)

                    fasta_number += 1
                    break

    return


if __name__ == "__main__":
    args = parse_args()

    # Establish globals if provided
    REGION_SIZE = args.gap
    UPSTREAM_ANCHOR = args.upstream
    DOWNSTREAM_ANCHOR = args.downstream

    main(args)
