"""
    Annotation to BED | RPINerd, 12/09/24

    Given the input from:
        https://hgdownload.soe.ucsc.edu/goldenPath/hgXX/database/rmsk.txt.gz,
    convert the file to a simple BED format with a summary description 4th column

    Input file columns:
        bin|swScore|milliDiv|milliDel|milliIns|GenoName|genoStart|genoEnd|genoLeft|strand|repName|repClass|repFamily|repStart|repEnd|repLeft|id

    BED info column:
        rep_type|total_len|unit_len|rep_len|rep_base
"""

import sys
from pathlib import Path

with Path.open(sys.argv[1], "r") as ref_file, Path.open("out.bed", "w") as output:
    for line in ref_file:
        columns = line.split("\t")

        genome_name = columns[5]
        genome_start = int(columns[6])
        genome_end = int(columns[7])
        rep_base = str(columns[10])
        rep_type = columns[11]
        total_len = genome_end - genome_start

        if rep_type == "Simple_repeat":
            unit_len = len(rep_base.strip("()n"))
            rep_len = int(total_len / unit_len)
        elif rep_type in {"Low_complexity", "Satellite"}:
            unit_len = "."
            rep_len = "."
        else:
            continue

        info_column = "|".join([rep_type, str(total_len), str(unit_len), str(rep_len), rep_base])
        bed_entry = "\t".join([genome_name, str(genome_start), str(genome_end), info_column])

        output.write(bed_entry + "\n")
