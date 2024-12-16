import sys

from Bio import SeqIO

r1_reads = {}
for record in SeqIO.parse(sys.argv[1], "fastq"):
    idr1 = record.id.split(" ")
    r1_reads[idr1[0]] = record

r2_reads = {}
for record in SeqIO.parse(sys.argv[2], "fastq"):
    idr2 = record.id.split(" ")
    r2_reads[idr2[0]] = record

r1_sync_file = open("R1_Sync.fastq", "w")
r2_sync_file = open("R2_Sync.fastq", "w")

for r2 in r2_reads.keys():
    if r1_reads.get(r2):
        SeqIO.write(r1_reads[r2], r1_sync_file, "fastq")
        SeqIO.write(r2_reads[r2], r2_sync_file, "fastq")
