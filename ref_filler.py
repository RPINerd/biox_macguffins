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

from simples import (look_backward_match, look_backward_miss,
                     look_forward_match, look_forward_miss)

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

    # Unpack refseq and set up some variables
    sequence = ref_seq
    gaps = []
    gap_start = None
    gap_end = None
    nt_total = len(sequence)

    # -print(id, sequence)
    idx = 0
    while idx < nt_total:

        if sequence[idx] != "-":  # Nothing to see here
            idx += REGION_SIZE - 1

        # A gap base is found at the current index
        else:

            # Find the start and end of the gap
            gap_start = look_backward_miss(sequence, idx, "-") + 1
            gap_end = look_forward_miss(sequence, idx, "-") - 1
            
            # Verify gap length is valid
            gap_len = gap_end - gap_start
            if gap_len >= REGION_SIZE:
                gaps.append((gap_start, gap_end))

            # Advance the index to the next base of the seq following the gap
            idx = gap_end + 1

    return gaps


def truncate_gaps(seq: tuple[str]) -> tuple[str]:
    """
    Truncate a given sequence to remove all gaps from both ends. 

    Such that:
    ATG-TTACTA-GTAGATAGCTGTCGTAGCGATGCAGTCGGCG-GCGAT-GGA--GATGCGA
    Would become:
    GTAGATAGCTGTCGTAGCGATGCAGTCGGCG

    :param tuple[str]: seq: Sequence to truncate
    :rtype: str - Truncated sequence
    """

    # Look backwards from the point where the upstream anchor joins to find the first gap
    upstream_first_gap = look_backward_match(seq, UPSTREAM_ANCHOR, "-") + 1

    # Look forward from the point where the downstream anchor joins to find the first gap
    downstream_first_gap = look_forward_match(seq, len(seq) - DOWNSTREAM_ANCHOR, "-") - 1

    #! Debug
    print(f"Upstream gap: {upstream_first_gap}, Downstream gap: {downstream_first_gap}")
    print("".join(seq[upstream_first_gap:downstream_first_gap]))

    return "".join(seq[upstream_first_gap:downstream_first_gap])


def generate_patch(tracks: dict[str, tuple[str]], start: int, end: int, output_dir: str, fasta_number: int) -> bool:

    start_upstream = start - UPSTREAM_ANCHOR
    end_downstream = end + DOWNSTREAM_ANCHOR

    for id, seq in tracks.items():
        if all(nt != "-" for nt in seq[start:end]):

            # TODO figure out how to best name these subsequences
            # Generate a fasta entry for the expansion
            fasta_header = f">ExpansionSEQ {fasta_number}"

            # TODO try except for when the downstream/upstream are outside of the range, default to just the start or end

            #! Debug
            print(f"Filler from: {id.strip()} with seq {"".join(seq[start:end])}")

            # If gaps are present in the upstream or downstream regions, truncate the sequence
            upstream_gaps = any(nt == "-" for nt in seq[start_upstream:start])
            downstream_gaps = any(nt == "-" for nt in seq[end:end_downstream])
            if upstream_gaps or downstream_gaps:
                print("Truncating gaps")
                fasta_seq = truncate_gaps(seq[start_upstream:end_downstream])
            else:
                fasta_seq = "".join(seq[start_upstream:end_downstream])

            #! Debug
            print(f"Seq including anchors:\n{fasta_seq}")

            # Write the fasta entry to a file, if the file already exists, overwrite it
            with open(os.path.join(output_dir, f"Expansion_{fasta_number}.fasta"), "w") as patch:
                patch.write(f"{fasta_header}\n{fasta_seq}\n")

            return True
    
    return False


def main(args) -> None:

    fasta_file = args.fasta
    output_directory = args.output

    # Verify that the output directory exists, create it if it does not
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
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
    fasta_number = 1
    for region in missing_regions:
        start, end = region
        if generate_patch(tracks, start, end, output_directory, fasta_number):
            fasta_number += 1

    return


if __name__ == "__main__":
    args = parse_args()

    # Establish globals if provided
    REGION_SIZE = args.gap
    UPSTREAM_ANCHOR = args.upstream
    DOWNSTREAM_ANCHOR = args.downstream

    main(args)
