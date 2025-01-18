"""
    CutAdapt on Directory | RPINerd, 01/18/25

    Tiny wrapper to run a cutadapt command, only for a very specific instance at this point
"""

import os
import re
from pathlib import Path

from configs import CUTADAPT_DIR

files = Path.iterdir(CUTADAPT_DIR)
adapters = CUTADAPT_DIR + "mmadapt.fasta"

for file in files:
    # TODO set to be universal
    file_info = re.match(r"(MMV2-R[0-9]{1,2}_S[0-9]{1,2}_L00[1234]_R)([12])(_001.fastq).gz", file)

    if not file_info:
        continue

    sample = file_info.group(1)
    read = file_info.group(2)
    suffix = file_info.group(3)

    if read == "2":
        continue

    print(file)
    print(sample + " | " + read + " | " + suffix + "\n")

    paired_read = sample + "2" + suffix + ".gz"
    short_r1 = "short_" + sample + str(read) + suffix
    short_r2 = "short_" + sample + "2" + suffix
    trimmed_r1 = "trimmed_" + sample + str(read) + suffix
    trimmed_r2 = "trimmed_" + sample + "2" + suffix

    header = 'echo "' + sample + '" >> reports.tsv'
    os.system(header)
    cut_cmd = f"cutadapt -j 12 --report=minimal -n 3 -m 25 -a file:{adapters} -A file:{adapters}    \
        --too-short-output {short_r1} --too-short-paired-output {short_r2} --pair-filter=first      \
        -o {trimmed_r1} -p {trimmed_r2} {file} {paired_read} >> reports.tsv"

    print("cutadapt:\n" + cut_cmd + "\n")

    os.system(cut_cmd)
