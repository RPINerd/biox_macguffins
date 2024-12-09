"""
"""
import sys
from pathlib import Path

import requests

server = "https://rest.ensembl.org"

with Path.open(sys.argv[1]) as id_list:
    for line in id_list:
        enst = line.strip().split(".")[0]
        ext = f"/sequence/id/{enst}?type=cdna"
        r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        print(r.text)
