"""
    cd-hit Cluster Splitter | RPINerd, 2021

    A cleanup script to break down output from cd-hit program
"""

import re
import sys

clstr_file = open(sys.argv[0])
cur_clst = ""
clusters = {}
first = True
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

for key in clusters:
    if len(clusters[key]) == 3:
        triples.extend(clusters[key])
    elif len(clusters[key]) <= 2:
        remainder.extend(clusters[key])
    else:
        clstr_out_name = str(key).replace(" ", "_").strip() + ".txt"
        clstr_out = open(clstr_out_name, "w")
        for hit in clusters[key]:
            clstr_out.write(str(hit) + "\n")
        clstr_out.close()

t_file = open("triples.txt", "w")
for hit in triples:
    t_file.write(str(hit) + "\n")
t_file.close()

r_file = open("remainder.txt", "w")
for hit in remainder:
    r_file.write(str(hit) + "\n")
r_file.close()
