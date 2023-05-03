import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo

save each fasta into a struct

for each fasta in struct
    
    for each other fasta in struct

        alignment = AlignIO.read(sys.argv[1], 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        summary_align.dumb_consensus(float(sys.argv[2]))
        
        if align better
            save to best alignment
        else next fasta

    save best match to new struct

    