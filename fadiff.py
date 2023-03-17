import sys
from Bio import SeqIO

'''
    fadiff.py
    Chris Howard | 2021

    Tiny wrapper around SeqIO to generate the unique sequences from two fasta files
    Messy and slow..
'''

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

# matching_records = {}
# only_old = {}
# for existing_record in curr_db:
#     hit = False
#     for new_record in new_db:
#         if existing_record.seq.upper() == new_record.seq.upper():
#             matching_records[new_record.id] = "\t".join([existing_record.id, str(existing_record.seq), new_record.id, str(new_record.seq)])
#             hit = True
#             break
#     if not hit:
#         only_old[existing_record.id] = str(existing_record.seq)

# only_new = {}
# for record in new_db:
#     if record.id not in matching_records:
#         only_new[record.id] = str(record.seq)

# print("Only Existing")
# for item in only_old.items():
#     print("\t".join([item[0], str(item[1])]))

# print("Only New")
# for item in only_new.items():
#     print("\t".join([item[0], str(item[1])]))

# print("Matching")
# for item in matching_records.items():
#     print("\t".join([item[0], str(item[1])]))