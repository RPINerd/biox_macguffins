"""
    Validate Bisulfite | RPINerd, 08/23/23

    Simple script to double-check all the primers desigend by the Bisulfite pipeline
    and make sure none create tabix hits (i.e. contain an SNP with frequency > 0.01)

"""

import os
import sys

from classes import Primer
from internalconfigs import TABIX_PATH

# TODO Merge this with crispr validate, basically identical!
assayfile = open(sys.argv[1], "r")
primers = []
pnum = 1

print("Reading input file...")
firstline = True
for line in assayfile:
    data = line.split(";")

    if firstline:
        firstline = False
        continue

    if data[23].strip():
        print("data23 " + data[23])
        continue

    p = Primer(data[2], pnum, data[4], data[3])
    primers.append(p)
    pnum += 1
print(f"Input file parsed! Found {len(primers)} primers.")

print(f"Using {TABIX_PATH} for index file.")
snp_found = 0
for primer in primers:
    tabix_query = f"tabix {TABIX_PATH} {primer.chr}:{primer.start}-{primer.end}"
    stream = os.popen(tabix_query)
    results = stream.readlines()
    for i in results:
        result_tabs = i.split("\t")
        if float(result_tabs[3]) <= 0.01:
            continue
        else:
            snp_found += 1
            print(f"SNP found on primer {primer.number}!")
            print(f"Primer loc: {primer.start}-{primer.end}  |  SNP at {result_tabs[1]}/{result_tabs[2]}\n")

print("All primers checked...")
print(f"{snp_found} SNP's found!")
