
file = open("base_content.tsv", "r")
clean = open("seq_list_cleaned.tsv", "w")
for line in file:

    if line.startswith('#'):
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
    elif ( a >= 75.00 ) | (t >= 75.00) | (c >= 75.00) | (g >= 75.00) | \
         (at >= 86.00) | (ac >= 86.00) | (ag >= 86.00) | (tc >= 86) | (tg >= 86) | (cg >= 86):
        continue
    else:
        clean.write(line)

file.close()
clean.close()