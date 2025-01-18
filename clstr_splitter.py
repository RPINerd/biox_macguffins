"""
    cd-hit Cluster Splitter | RPINerd, 01/18/25

    A cleanup script to break down output from cd-hit program
"""

import re
import sys
from pathlib import Path

cur_clst = ""
clusters = {}
first = True
with Path.open(sys.argv[0], "r") as clstr_file:
    for line in clstr_file:
        if first:
            first = False
            cur_clst = line.strip(">")
            clusters[cur_clst] = []
            continue

        if line.startswith(">"):
            cur_clst = line.strip(">")
            clusters[cur_clst] = []

        else:
            reg = r">(chr.+:[0-9]+-[0-9]+)\.\."
            clusters[cur_clst].append(re.search(reg, line).group(1))


print("Size of clusters: " + str(len(clusters)) + "\n")

triples = []
remainder = []

for key in clusters.items():
    if len(clusters[key]) == 3:
        triples.extend(clusters[key])
    elif len(clusters[key]) <= 2:
        remainder.extend(clusters[key])
    else:
        clstr_out_name = str(key).replace(" ", "_").strip() + ".txt"
        with Path.open(clstr_out_name, "w") as clstr_out:
            for hit in clusters[key]:
                clstr_out.write(str(hit) + "\n")

with Path.open("triples.txt", "w") as t_file:
    for hit in triples:
        t_file.write(str(hit) + "\n")

with Path.open("remainder.txt", "w") as r_file:
    for hit in remainder:
        r_file.write(str(hit) + "\n")
