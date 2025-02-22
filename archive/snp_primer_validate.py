"""
    Primer SNP Validation | RPINerd, 08/23/23

    Simple script to double-check primers designed by either the Bisulfite or Crispr
    pipeline to make sure none create tabix hits (i.e. contain an SNP with frequency > 0.01)

    Usage: python [assay file] [b/c]
"""

import os
import sys
from pathlib import Path

from internal.internalconfigs import TABIX_PATH

from classes import Primer

SNP_THRESHOLD = 0.01

primers: list[Primer] = []
pnum: int = 1

print("Reading input file...")
with Path.open(sys.argv[1]) as assayfile:
    firstline = True
    for line in assayfile:
        data = line.split(";")

        if firstline:
            firstline = False
            continue

        # Bisulfite primers input
        if sys.argv[2] == "b":
            if data[23].strip():
                print("data23 " + data[23])
                continue
            p = Primer(data[2], pnum, data[4], data[3])
            primers.append(p)
            pnum += 1

        # Crispr primers input
        elif sys.argv[2] == "c":
            fp = Primer(data[0], pnum, data[4], data[5])
            rp = Primer(data[0], pnum, data[11], data[10])
            primers.append(fp)
            primers.append(rp)
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
        if float(result_tabs[3]) <= SNP_THRESHOLD:
            continue
        else:
            snp_found += 1
            print(f"SNP found on primer {primer.number}!")
            print(f"Primer loc: {primer.start}-{primer.end}  |  SNP at {result_tabs[1]}/{result_tabs[2]}\n")

print("All primers checked...")
print(f"{snp_found} SNP's found!")
