import sys

'''
    ann_to_bed.py
    Chris Howard | 04/04/23
    
    Given the input from https://hgdownload.soe.ucsc.edu/goldenPath/hgXX/database/rmsk.txt.gz, convert the file to a simple BED format with a summary description 4th column

    Input file columns:
        bin | swScore | milliDiv | milliDel | milliIns | GenoName	genoStart | genoEnd | genoLeft | strand | repName | repClass	repFamily | repStart | repEnd | repLeft | id

    BED info column:
        repType|totalLen|unitLen|repLen|repBase

    
'''

try:
    ref_file = open(sys.argv[1], "r")
except FileNotFoundError:
    print("Reference file not found!")
else:
    with open("out.bed", "w") as output:
        for line in ref_file:
            columns = ref_file.split("\t")

            repType     = columns[11]
            if repType not in ["Simple_repeat", "Low_complexity", "Satellite"]:
                continue
            genomeName  = columns[5]
            genomeStart = int(columns[6])
            genomeEnd   = int(columns[7])
            totalLen    = genomeEnd - genomeStart
            repBase     = str(columns[10])
            unitLen     = len(repBase.strip("()n"))
            repLen      = int(totalLen / unitLen)

            info_column = "|".join([repType, totalLen, unitLen, repLen, repBase])
            bed_entry = "\t".join([genomeName, genomeStart, genomeEnd, info_column])

            output.write(bed_entry + "\n")
