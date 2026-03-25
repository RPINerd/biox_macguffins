"""Example configs file to set up paths and variables"""

from pathlib import Path

from Bio import Entrez

Entrez.email = "rpimule@gmail.com"

REFSEQ_CACHE_DIR = Path(__file__).parent.parent / "refseqs"

TABIX_PATH = "/home/user/tabix/"
MIN_LENGTH_COMPLEXITY = 50
