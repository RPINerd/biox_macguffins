"""
    Fasta Diff | RPINerd, 2021

    Tiny wrapper around SeqIO to generate the unique sequences from two fasta files
    Messy and slow..
"""

import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

curr_db: list[SeqRecord] = list(SeqIO.parse(sys.argv[1], "fasta"))
new_db: list[SeqRecord] = list(SeqIO.parse(sys.argv[2], "fasta"))

uniq_seqs = {}
for record in curr_db:
    if str(record.seq).lower() in uniq_seqs:
        uniq_seqs[str(record.seq).lower()].append(str(record.id))
    else:
        uniq_seqs[str(record.seq).lower()] = [str(record.id)]

for seq, old_entries in uniq_seqs.items():
    uniq_seqs[seq] = [",".join(old_entries)]

for record in new_db:
    if str(record.seq).lower() in uniq_seqs:
        uniq_seqs[str(record.seq).lower()].append(str(record.id))
    else:
        uniq_seqs[str(record.seq).lower()] = ["N/A", str(record.id)]

for seq, ids in uniq_seqs.items():
    print("{}\t{}".format(seq, "\t".join(ids)))
