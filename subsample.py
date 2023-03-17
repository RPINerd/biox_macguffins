#!/usr/bin/python3
import sys
import random
import os
import argparse
import subprocess
import re
from pathlib import Path

'''
    subsample.py
    Chris Howard | 05/27/21

    Wrapper around SeqTK's sample function to generate a subset fastq file based on the file given by the user and the number of subsamples requested.
'''

# TODO allow parsing of a file to bulk subsample
def parse_tsv():
  return


def main():

    ### Argument parsing
    parser = argparse.ArgumentParser(description="Generate a fastq file with a random subset of reads from the given input file. Default subset is 2000.")
    parser.add_argument("-f", "--file", help="File, or tsv file containing reads you wish to subset", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    #parser.add_argument("-o", "--output", help="Designates a user-defined output file")
    parser.add_argument("-n", "--number", help="Number of reads desired in the subset", type=int, default=2000)
    args = parser.parse_args()


    ### Processing of input file
    # Make sure it exists
    if not os.path.isfile(args.file):
        print('Input file does not exist!')
        exit
    else:
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

    ### Prepare output files
    outfile_R1 = sample_id + "_R1." + str(args.number) + ".fastq"
    outfile_R2 = sample_id + "_R2." + str(args.number) + ".fastq"
    if args.verbose:
        print("Writing out to files: " + outfile_R1 + ", " + outfile_R2)

    seed = random.randint(1,999)
    if args.verbose:
        print("Seed value for this run is: " + str(seed))

    # Execute the call
    # TODO does seqtk handle/alert when subsample is >= input read amount?
    sub1 = "seqtk sample -s" + str(seed) + " " + read1_file + " " + str(args.number) + " > " + outfile_R1
    if args.verbose:
        print("Read 1 call generated and executing: \n" + sub1)
    os.system(sub1)
    sub2 = "seqtk sample -s" + str(seed) + " " + read2_file + " " + str(args.number) + " > " + outfile_R2
    if args.verbose:
        print("Read 2 call generated and executing: \n" + sub2)
    os.system(sub2)


if __name__ == "__main__":
    main()
