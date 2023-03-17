import os
import argparse
from posix import EX_OSFILE
import re
from pathlib import Path

'''
    grep_reads
    Chris Howard | 06/14/21

    Given a regex (or file of expressions) to look for, grep_reads will find all reads that contain said pattern in the R1 fastq files provided. Then it will write these desired reads to a new fastq file and then extract the matching paired reads from R2.
'''

def main():

    # Argument Parsing
    parser = argparse.ArgumentParser(description="Generate a paired fastq file with a subset of reads matching a user-provided regex.")
    parser.add_argument("-f", "--file", help="R1 fastq file to extract reads from", required=True)
    parser.add_argument("-r", "--regex", help="Regular expression of desired pattern contained in read", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    #parser.add_argument("-o", "--output", help="Designates a user-defined output file")
    args = parser.parse_args()


    ##Processing Inputs

    # Make sure input file exists
    input_file = args.file
    if not os.path.isfile(input_file):
        print('Input file does not exist!')
        exit

    # Check to see if it is a singular expression or a file
    regex = args.regex
    regex_file = False
    if os.path.isfile(regex):
        # Regex argument points to a file
        regex_file = True
        # patterns = []
        # for line in open(regex).read().splitlines():
        #     patterns.append(line)
        # regex = None

    # Gather a readset identifier from the provided R1 file
    sample_match = re.match("(.*)_R1(.*).fastq", input_file)
    id = sample_match.group(1)
    suffix = sample_match.group(2)
    if args.verbose:
        print( "Read ID: " + id )
        print( "Suffix: " + suffix )

    ## Collection and pairing of reads
    # Use provided pattern to collect R1 reads
    r1_out = id + "_R1_Matches.fastq"
    if regex_file:
        cmd = 'grep -A2 -B1 -f ' + regex + ' ' + input_file + ' > ' + r1_out
    else:
        cmd = 'grep -P -A2 -B1 "' + regex + '" ' + input_file + ' > ' + r1_out
    if args.verbose:
        print( "Finding read matches:\n" + cmd )
    os.system(cmd)

    # Clean up some grep artifacts
    cmd = "sed -i '/^--$/d' " + r1_out
    if args.verbose:
        print( "Cleaning matches (" + cmd + ")" )
    os.system(cmd)

    # Find matching R2 reads
    r2_input = id + "_R2" + suffix + ".fastq"
    r2_out = id + "_R2_Matches.fastq"
    cmd = 'for i in $(grep -Po "^@.*\s" ' + r1_out + '); do grep -A3 "$i" ' + r2_input + ' >> ' + r2_out + '; done;'
    if args.verbose:
        print( "Collecting paired ends:\n" + cmd )
    os.system(cmd)


if __name__ == "__main__":
    main()
