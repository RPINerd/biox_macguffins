"""
    Fasta Diff | RPINerd, 2021

    Tiny wrapper around SeqIO to generate the unique sequences from two fasta files
    Messy and slow..
"""

import sys

from Bio import SeqIO

curr_db = list(SeqIO.parse(sys.argv[1], "fasta"))
new_db = list(SeqIO.parse(sys.argv[2], "fasta"))

uniq_seqs = {}
for record in curr_db:
    if str(record.seq).lower() in uniq_seqs:
        uniq_seqs[str(record.seq).lower()].append(str(record.id))
    else:
        uniq_seqs[str(record.seq).lower()] = [str(record.id)]

for seq in uniq_seqs:
    old_entries = ",".join(uniq_seqs[seq])
    uniq_seqs[seq] = [old_entries]

for record in new_db:
    if str(record.seq).lower() in uniq_seqs:
        uniq_seqs[str(record.seq).lower()].append(str(record.id))
    else:
        uniq_seqs[str(record.seq).lower()] = ["N/A", str(record.id)]

for seq, ids in uniq_seqs.items():
    print("{}\t{}".format(seq, "\t".join(ids)))
