"""
    Pairwise Align | RPINerd, 08/22/23

    Given two input files in fasta format, align each pair of sequences and output the alignments into a single file

    Usage: python align.py <input_file_1> <input_file_2> <output_file>
    Example: python align.py input1.fasta input2.fasta output.fasta

    Note: This script requires the Biopython module: http://biopython.org/wiki/Download
  """

import sys

from Bio import Align, Seq, SeqIO

# Read in the input files
input_file_1 = sys.argv[1]
input_file_2 = sys.argv[2]
output_file = sys.argv[3]

# Create lists for the sequences and IDs in the input files
sequences_1 = []
sequences_2 = []
names_1 = []
names_2 = []
for record in SeqIO.parse(input_file_1, "fasta"):
    if len(record.seq) <= 30:
        seq = Seq.reverse_complement(record.seq)
    else:
        seq = record.seq
    sequences_1.append(seq)
    names_1.append(record.id)
for record in SeqIO.parse(input_file_2, "fasta"):
    sequences_2.append(record.seq)
    names_2.append(record.id)

# Set up an aligner
aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.open_gap_score = -2.0
aligner.extend_gap_score = -0.5
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0

# Create a list of the alignments
alignments = []
for i in range(len(sequences_1)):
    for j in range(len(sequences_2)):
        alignments.append(
            [
                f">{names_1[i]}::{names_2[j]}",
                aligner.align(sequences_1[i], sequences_2[j]),
            ]
        )

# Write the alignments to the output file
with open(output_file, "w") as output:
    for i in range(len(alignments)):
        output.write(alignments[i][0] + "\n")
        output.write(str(alignments[i][1][0][0]) + "\n")
        output.write(str(alignments[i][1][0][1]) + "\n")
        output.write("\n")
