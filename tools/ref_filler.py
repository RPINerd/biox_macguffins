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
from pathlib import Path

from macguffin_utils import look_backward_match, look_backward_miss, look_forward_match, look_forward_miss

# Minimum gap size in refseq to consider filling
REGION_SIZE = 150
# Bases up/down stream of patch to include
UPSTREAM_ANCHOR = 100
DOWNSTREAM_ANCHOR = 100


def parse_args() -> argparse.ArgumentParser:
    """Simple argument parser for the script"""
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

    return parser.parse_args()


# TODO make sure there is more than 1 sequence
# TODO input must be an alignment, i.e. have '-' characters
def validate_falign(file) -> str:
    """
    Validate the input fasta alignment file
    """
    message = ""

    return message


def parse_fasta(file: str) -> dict[str, tuple[str]]:
    """
    Parses a fasta alignment file into a dictionary of sequences.

    Args:
        file (str): Path to the fasta alignment file.

    Returns:
        dict[str, tuple[str]]: Dictionary of sequences with the ID as the key.
    """
    tracks = {}
    with Path.open(file, "r") as alignment:
        lines = alignment.readlines()

        ref = True
        for id, seq in zip(lines[::2], lines[1::2]):
            if ref:
                tracks["refseq"] = tuple(seq)
                ref = False
            else:
                tracks[id] = tuple(seq)
                # ! Debug - create a file to add all the insertions associated with the ID to in alignment format
                # with open(f"{id.strip('>').strip()}_insertions.fasta", "w") as insertions:
                #     insertions.write(f"{id.strip()}\n{seq}")

    return tracks


def parse_refseq(ref_seq: tuple[str]) -> list[tuple[int, int]]:
    """
    Parse the reference sequence to identify missing regions

    Args:
        ref_seq (tuple[str]): Tuple of the reference sequence

    Returns:
        list[tuple[int, int]]: List of tuples containing the start and end of missing regions
    """
    # Unpack refseq and set up some variables
    sequence = ref_seq
    gaps = []
    gap_start = None
    gap_end = None
    nt_total = len(sequence)

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


def generate_patch(tracks: dict[str, tuple[str]], start: int, end: int, output_dir: str, fasta_number: int) -> bool:
    """
    Generate a fasta file for the missing region

    Args:
        tracks (dict[str, tuple[str]]): Dictionary of sequences
        start (int): Start of the missing region
        end (int): End of the missing region
        output_dir (str): Directory to write the fasta files to
        fasta_number (int): Number to append to the fasta file name

    Returns:
        bool: True if a patch was generated, False otherwise
    """
    upstream_end = start - 1
    upstream_start = upstream_end - UPSTREAM_ANCHOR
    downstream_start = end + 1
    downstream_end = downstream_start + DOWNSTREAM_ANCHOR

    for id, seq in tracks.items():
        if any(nt == "-" for nt in seq[start:end]):
            continue

        filler_seq = "".join(seq[start:end])

        # TODO account for when the downstream/upstream are outside of the range, default to just the start or end
        upstream_seq = "".join(seq[upstream_start:upstream_end])
        downstream_seq = "".join(seq[downstream_start:downstream_end])

        # If gaps are present in the upstream or downstream regions, truncate the sequence
        if "-" in upstream_seq:
            upstream_cut_point = look_backward_match(upstream_seq, len(upstream_seq) - 1, "-")
            upstream_seq = upstream_seq[upstream_cut_point:]
        if "-" in downstream_seq:
            downstream_cut_point = look_forward_match(downstream_seq, 0, "-") + 1
            downstream_seq = downstream_seq[:downstream_cut_point]

        # Extract the sequence to be written to the fasta file
        fasta_filename = f"{fasta_number}_{id.strip(">").strip()}_{upstream_start}-{downstream_end}.fasta"
        fasta_header = f"{id.strip()}:{upstream_start}-{downstream_end}"
        fasta_seq = "".join([upstream_seq, filler_seq, downstream_seq])

        # Write the fasta entry to a file, if the file already exists, overwrite it
        with Path.open(Path.joinpath(output_dir, fasta_filename), "w") as patch:
            patch.write(f"{fasta_header}\n{fasta_seq}\n")

        return True

    return False


def main(args: argparse.Namespace) -> None:
    """
    Driver function for the script

    Args:
        args (argparse.Namespace): Arguments passed to the script

    Returns:
        None

    """
    fasta_file = args.fasta
    output_directory = args.output

    # Verify that the output directory exists, create it if it does not
    if not Path.exists(output_directory):
        Path.mkdir(output_directory, parents=True)

    # TODO validate fasta file instead of assuming valid input
    invalid = validate_falign(fasta_file)
    if invalid:
        raise f"Invalid fasta alignment input: {invalid}"

    # Extract the individual sequences as sets
    tracks = parse_fasta(fasta_file)

    # Process the reference sequence to establish relevant gaps
    ref_seq = tracks.pop("refseq")
    missing_regions = parse_refseq(ref_seq)

    # Process the other sequences to extract the missing regions
    fasta_number = 1
    for region in missing_regions:
        start, end = region
        if generate_patch(tracks, start, end, output_directory, fasta_number):
            fasta_number += 1


if __name__ == "__main__":
    args = parse_args()

    # Overwrite globals if provided
    REGION_SIZE = args.gap
    UPSTREAM_ANCHOR = args.upstream
    DOWNSTREAM_ANCHOR = args.downstream

    main(args)
