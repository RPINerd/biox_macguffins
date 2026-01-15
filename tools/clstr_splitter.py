"""
    cd-hit Cluster Splitter | RPINerd, 01/18/25

    A cleanup script to break down output from cd-hit program
"""

import re
import sys
from pathlib import Path

cur_clst = ""
clusters: dict[str, list[str]] = {}
first = True
clst_file = Path(sys.argv[0])
if not clst_file.exists():
    raise FileNotFoundError(f"Cluster file {clst_file} does not exist!")

with Path.open(clst_file) as infile:
    for line in infile:
        if first:
            first = False
            cur_clst = line.strip(">")
            clusters[cur_clst] = []
            continue

        if line.startswith(">"):
            cur_clst = line.strip(">")
            clusters[cur_clst] = []

        else:
            regex = r">(chr.+:[0-9]+-[0-9]+)\.\."
            hit = re.search(regex, line)
            if not hit:
                print(f"Line {line} does not appear to have expected structure")
                continue
            region = hit.group(1)
            clusters[cur_clst].append(region)


print("Size of clusters: " + str(len(clusters)) + "\n")

triples: list = []
remainder: list = []
triples_file: Path = Path("triples.txt")
remainder_file: Path = Path("remainder.txt")

for key in clusters.items():
    if len(clusters[key]) == 3:
        triples.extend(clusters[key])
    elif len(clusters[key]) <= 2:
        remainder.extend(clusters[key])
    else:
        clstr_out_name = Path(clst_file).parent / (str(key).replace(" ", "_").strip() + ".txt")
        with Path.open(clstr_out_name, "w") as clstr_out:
            for hit in clusters[key]:
                clstr_out.write(str(hit) + "\n")

with triples_file.open("w", encoding="utf-8") as t_file:
    for hit in triples:
        t_file.write(str(hit) + "\n")

with remainder_file.open("w", encoding="utf-8") as r_file:
    for hit in remainder:
        r_file.write(str(hit) + "\n")
