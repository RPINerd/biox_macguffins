"""
    Window Shopper | RPINerd, 06/14/23

    Given a fasta file, step through all the sequences and break down each one into a collection of X-base windows.
"""

import argparse

from Bio import SeqIO

import assembly


def main(args) -> None:
    window_size = int(args.window)
    window_dict = {}
    for record in SeqIO.parse(args.fasta, "fasta"):
        seq = record.seq
        step = int(args.step)
        start = 0
        while start < len(seq) - window_size:
            end = start + window_size
            window = seq[start:end]
            try:
                window_dict[window] += 1
            except KeyError:
                window_dict[window] = 1
            start += step

    if args.output:
        i = 0
        with open(args.output, "w") as fa_file:
            for seq, count in window_dict.items():
                fa_file.write(f">{i}-{count}\n")
                fa_file.write(f"{seq}\n")
                i += 1
    else:
        i = 0
        for seq, count in window_dict.items():
            print(f">{i}-{count}")
            print(f"{seq}")
            i += 1


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--fasta", help="Input fasta file", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    parser.add_argument("-o", "--output", help="Designates a user-defined output file. (Default prints to stdout)")
    parser.add_argument("-w", "--window", help="Size of the sliding window (Default = 30bp)", default=30)
    parser.add_argument("-a", "--assembly", help="Perform assembly instead of window slicing.", action="store_true")
    parser.add_argument("-s", "--step", help="How many steps to take per window iteration", default=1)
    args = parser.parse_args()

    if args.assembly:
        s = assembly.SequenceAssembler()

        # read data
        with open(args.fasta) as data:
            s.read_fasta(data)

        # run assembly
        s.assemble()

        # access path and resulting sequence

        print(s.path)
        print(s.sequence)

    else:
        main(args)
