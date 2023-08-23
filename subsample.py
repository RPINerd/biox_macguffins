"""
    SeqTK Subsampling | RPINerd, 08/23/23

    Wrapper around SeqTK's sample function to generate a subset fastq file based
    on the file given by the user and the number of subsamples requested.
"""

import argparse
import os
import random
import re


# TODO allow parsing of a file to bulk subsample
def parse_tsv():
    return


def arg_parse():
    parser = argparse.ArgumentParser(
        description="Generate a fastq file with a random subset of reads from the given input file"
    )
    parser.add_argument("-f", "--file", help="File, or tsv file containing reads you wish to subset", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    # parser.add_argument("-o", "--output", help="Designates a user-defined output file")
    parser.add_argument(
        "-n", "--number", help="Number of reads desired in the subset (default: 2000)", type=int, default=2000
    )
    args = parser.parse_args()

    return args


def main(args):
    ## Processing of input file
    # Make sure it exists
    assert os.path.isfile(args.file), "Input file does not exist!"

    if args.verbose:
        print("Input file name: " + args.file)

    # User provided *.tsv
    if re.search(r".tsv$", args.file):
        # Process the list of desired reads
        read_files = parse_tsv

    # User provided single fastq
    elif re.search(r"fastq.gz$", args.file):
        file_reg = re.search(r"(^.*)_R[12](.*fastq.gz)", args.file)
        sample_id = file_reg.group(1)
        suffix = file_reg.group(2)
        read1_file = sample_id + "_R1" + suffix
        read2_file = sample_id + "_R2" + suffix
        if args.verbose:
            print("Sample ID: " + sample_id)
            print("File Suffix: " + suffix)
            print("Read files: " + read1_file + ", " + read2_file)

    ## Prepare output files
    outfile_R1 = f"{sample_id}_R1.{args.number}.fastq"
    outfile_R2 = f"{sample_id}_R2.{args.number}.fastq"
    if args.verbose:
        print("Writing out to files: " + outfile_R1 + ", " + outfile_R2)

    seed = random.randint(1, 999)
    if args.verbose:
        print("Seed value for this run is: " + str(seed))

    # Execute the call
    # TODO does seqtk handle/alert when subsample is >= input read amount?
    sub1 = f"seqtk sample -s {seed} {read1_file} {args.number} > {outfile_R1}"
    if args.verbose:
        print("Read 1 call generated and executing: \n" + sub1)
    os.system(sub1)
    sub2 = f"seqtk sample -s {seed} {read2_file} {args.number} > {outfile_R2}"
    if args.verbose:
        print("Read 2 call generated and executing: \n" + sub2)
    os.system(sub2)


if __name__ == "__main__":
    args = arg_parse()
    main(args)
