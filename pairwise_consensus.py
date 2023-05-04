import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo

def main(file):

    records = []
    try:
        for record in SeqIO.parse(file, "fasta"):
            records.append(record)
    except Exception as e:
        print(f"Error on record parsing: {e}")
        exit

    for sequence in records:
        reference_list = records.remove(sequence)
        for ref_sequence in reference_list:

            alignment = AlignIO.align()
            summary_align = AlignInfo.SummaryInfo(alignment)
            summary_align.dumb_consensus(float(sys.argv[2]))
            
            if align better
                save to best alignment
            else next fasta

        save best match to new struct

if __name__ == "__main__":
    input_fasta = sys.argv[1]
    main(input_fasta)