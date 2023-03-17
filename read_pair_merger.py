import sys
import os
import re
from Bio import SeqIO

r1_file = sys.argv[1]
r2_file = sys.argv[2]
result_file = sys.argv[3]

#! Debug Print
#print(r1_file)
#print(r2_file)

reads_1 = []
reads_2 = []
r1_count = 0
r2_count = 0
t_run = 0

# Read 1 ingest
for record in SeqIO.parse(r1_file, "fastq"):
    reads_1.append(record.seq)
    r1_count += 1

# Read 2 ingest
for record in SeqIO.parse(r2_file, "fastq"):
    reads_2.append(record.seq)
    r2_count += 1

    # Read 2 contains a T run of 20 or more, count for final percentage
    if re.search("T{20,}", str(record.seq)):
        t_run += 1

if (r1_count != r2_count):

    print("Unequal count! Perhaps you mismatched R1/R2 files?")

else:

    # Open ouput csv for report
    out_file = open(result_file, "w")

    # Write data
    r = 0
    header = True
    while r < r1_count:

        # First line info
        if header:
            header = False
            out_file.write("Read 1,Read 2,%R2 >= 20T : " + str( (t_run / r2_count) * 100 ) + "\n")

            #! Debug Print
            #print("Read 1\tRead 2\t%R2 >= 20T : " + str( (t_run / r2_count) * 100 ))

            next

        #! Debug Print
        #print("R1: " + reads_1[r] + "\nR2: " + reads_2[r] + "-------------------")

        out_file.write( str(reads_1[r]) + "," + str(reads_2[r]) + "\n")

        r += 1