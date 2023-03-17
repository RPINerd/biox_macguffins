import sys
import os

# Internal
from configs import *
from classes import Primer

'''
    bisulf_validate
    Chris Howard | 03/17/23
    
    Simple script to double-check all the primers desigend by the Bisulfite pipeline and make sure none create tabix hits (i.e. contain an SNP with frequency > 0.01)
'''

assayfile = open(sys.argv[1],"r")
primers = []
pnum = 1

print("Reading input file...")
firstline = True
for line in assayfile:

    data = line.split(";")

    if ( firstline ):
        firstline = False
        continue

    if ( data[23].strip() ):
        print("data23 " + data[23])
        continue

    p = Primer(data[2], pnum, data[4], data[3])
    primers.append(p)
    pnum += 1
print("Input file parsed! Found " + str(len(primers)) + " primers.")

print("Using " + TABIX_PATH + " for index file.")
snp_found = 0
for primer in primers:
    tabix_query = "tabix " + TABIX_PATH + " " + primer.chr + ":" + primer.start + "-" + primer.end
    stream = os.popen(tabix_query)
    results = stream.readlines()
    for i in results:
        result_tabs = i.split('\t')
        if (float(result_tabs[3]) <= 0.01):
            continue
        else:
            snp_found += 1
            print("SNP found on primer " + str(primer.number) + "!")
            print("Primer loc: " + str(primer.start) + "-" + str(primer.end) + "  |  SNP at " + result_tabs[1] + "/" + result_tabs[2] + "\n")

print("All primers checked...")
if ( snp_found > 0 ):
    print( str(snp_found) + " SNP's found!")
else:
    print("No SNP's found!")
