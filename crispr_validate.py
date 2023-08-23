"""
    Validate Crispr | RPINerd, 03/17/23

    Simple script to double-check all the primers desigend by the CRISPR pipeline
    and make sure none create tabix hits (i.e. contain an SNP with frequency > 0.01)
"""

import os
import sys

from classes import Primer
from configs import TABIX_PATH

assayfile = open(sys.argv[1], "r")
primers = []
pnum = 1

print("Reading input file...")
for line in assayfile:
    data = line.split(";")

    if data[0] == "Chrom":
        continue

    fp = Primer(data[0], pnum, data[4], data[5])
    rp = Primer(data[0], pnum, data[11], data[10])
    primers.append(fp)
    primers.append(rp)
    pnum += 1
print("Input file parsed! Found " + str(len(primers)) + " primers.")

print("Using " + TABIX_PATH + " for index file.")
snp_found = 0
for primer in primers:
    tabix_query = "tabix " + TABIX_PATH + " " + primer.chr + ":" + primer.start + "-" + primer.end
    stream = os.popen(tabix_query)
    results = stream.readlines()
    for i in results:
        result_tabs = i.split("\t")
        if float(result_tabs[3]) <= 0.01:
            continue
        else:
            snp_found += 1
            print("SNP found on primer " + str(primer.number) + "!")
            print(f"Primer loc: {primer.start}-{primer.end}  |  SNP at {result_tabs[1]}/{result_tabs[2]}\n")

print("All primers checked...")
print(f"{snp_found} SNP's found!")
