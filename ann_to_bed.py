"""
    Annotation to BED | RPINerd, 08/23/23

    Given the input from:
        https://hgdownload.soe.ucsc.edu/goldenPath/hgXX/database/rmsk.txt.gz,
    convert the file to a simple BED format with a summary description 4th column

    Input file columns:
        bin|swScore|milliDiv|milliDel|milliIns|GenoName|genoStart|genoEnd|genoLeft|strand|repName|repClass|repFamily|repStart|repEnd|repLeft|id

    BED info column:
        repType|totalLen|unitLen|repLen|repBase
"""

import sys

with open(sys.argv[1], "r") as ref_file:
    with open("out.bed", "w") as output:
        for line in ref_file:
            columns = line.split("\t")

            genomeName = columns[5]
            genomeStart = int(columns[6])
            genomeEnd = int(columns[7])
            repBase = str(columns[10])
            repType = columns[11]
            totalLen = genomeEnd - genomeStart

            if repType == "Simple_repeat":
                unitLen = len(repBase.strip("()n"))
                repLen = int(totalLen / unitLen)
            elif repType in ["Low_complexity", "Satellite"]:
                unitLen = "."
                repLen = "."
            else:
                continue

            info_column = "|".join([repType, str(totalLen), str(unitLen), str(repLen), repBase])
            bed_entry = "\t".join([genomeName, str(genomeStart), str(genomeEnd), info_column])

            output.write(bed_entry + "\n")
