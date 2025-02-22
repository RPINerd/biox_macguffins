"""
    Fasta Unique | RPINerd, 2021

    Tiny wrapper around SeqIO to generate the unique sequences from a single fasta file
"""

import sys

from Bio import SeqIO

new_db = list(SeqIO.parse(sys.argv[1], "fasta"))

uniq_seqs = {}
for record in new_db:
    if record.seq in uniq_seqs:
        uniq_seqs[str(record.seq)].append(str(record.id))
    else:
        uniq_seqs[str(record.seq)] = [str(record.id)]

for entry in uniq_seqs:
    ids = "\t".join(uniq_seqs[entry])
    print("\t".join([entry, ids]))
