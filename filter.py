"""
    Base Content Filter | RPINerd, ??

    Input is a tab separated list of sequences with base content metrics:
    seq|?|gc-poly|a|t|c|g
    With each entry done in percent values

    Script will output a cleaned tsv file that removes all sequences not meeting 
    desired thresholds
"""

file = open("base_content.tsv", "r")
clean = open("seq_list_cleaned.tsv", "w")
for line in file:
    if line.startswith("#"):
        continue

    cols = line.split()

    gcp = float(cols[2])
    a = float(cols[3])
    t = float(cols[4])
    c = float(cols[5])
    g = float(cols[6])

    at = a + t
    ac = a + c
    ag = a + g
    tc = t + c
    tg = t + g
    cg = c + g

    if gcp >= 75.00:
        continue
    elif (
        (a >= 75.00)
        or (t >= 75.00)
        or (c >= 75.00)
        or (g >= 75.00)
        or (at >= 86.00)
        or (ac >= 86.00)
        or (ag >= 86.00)
        or (tc >= 86)
        or (tg >= 86)
        or (cg >= 86)
    ):
        continue
    else:
        clean.write(line)

file.close()
clean.close()
