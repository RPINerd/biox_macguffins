import argparse
import sys
import re

'''
    cigar_filter.py
    Chris Howard | 03/21/23
    
    For now just a barebones script to filter out reads based on a single CIGAR characteristic. Potential here to develop into something much more full featured.
'''


def main() -> None:
    
    out = open(sys.argv[2], "w")

    for line in open(sys.argv[1], "r"):
        if line.startswith("@"):
            continue
        cols = line.split("\t")
        # cols[5] is cigar
        filter = re.match("([0-9]+)M", cols[5])
        if filter and int(filter[1]) >= 50:
            print(cols[2], cols[3], cols[5])
            out.write("\t".join([cols[2], cols[3], cols[5]]) + "\n")

if __name__ == "__main__":
    main()