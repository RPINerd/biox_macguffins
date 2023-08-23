import re
import sys

from Bio import SeqIO

r1_file = sys.argv[1]
r2_file = sys.argv[2]
result_file = sys.argv[3]

reads_1 = []
reads_2 = []
t_run = 0

# Read 1 ingest
for record in SeqIO.parse(r1_file, "fastq"):
    reads_1.append(record.seq)

# Read 2 ingest
for record in SeqIO.parse(r2_file, "fastq"):
    reads_2.append(record.seq)

    # Read 2 contains a T run of 20 or more, count for final percentage
    if re.search("T{20,}", str(record.seq)):
        t_run += 1

r1_count = len(reads_1)
r2_count = len(reads_2)
if r1_count != r2_count:
    print("Unequal count! Perhaps you mismatched R1/R2 files?")

else:
    # Open ouput csv for report
    out_file = open(result_file, "w")

    # Write data
    polyT_percent = (t_run / r2_count) * 100
    out_file.write(f"Read 1,Read 2,%R2 >= 20T : {polyT_percent}\n")
    for r in range(r1_count):
        out_file.write(f"{reads_1[r]},{reads_2[r]}\n")

        r += 1
