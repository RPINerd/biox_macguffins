"""
    Acquisition Functions | RPINerd, 03/25/26

    Functions for acquiring sequences or other refernce data from public databases, such as NCBI, Ensembl, etc.
"""

import logging

import requests
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from biopython_plus import extract_exons
from configs import REFSEQ_CACHE_DIR

ENSEMBL_REST = "https://rest.ensembl.org"

logging.basicConfig(
    level=logging.DEBUG,
    format=("%(asctime)s | %(levelname)-7s | %(module)-10s | %(lineno)-4d | %(message)s"),
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


def ensembl_fetch_from_ids(id_list: list[str]) -> dict[str, str]:
    """
    Fetch sequences from Ensembl REST API given a list of transcript IDs

    ! Needs to be reworked, hits the rest API for each ID, not polite

    Args:
        id_list (list[str]): List of Ensembl transcript IDs (e.g. ENST00000367770)

    Returns:
        dict[str, str]: Dictionary mapping transcript IDs to their corresponding FASTA sequences
    """
    sequences = {}
    for line in id_list:
        enst = line.strip().split(".")[0]

        # Check cache first
        cache_file = REFSEQ_CACHE_DIR / f"{enst}.fa"
        if cache_file.exists():
            with cache_file.open("r") as f:
                sequences[enst] = f.read()
            continue

        ext = f"/sequence/id/{enst}?type=cdna"
        req = requests.get(ENSEMBL_REST + ext, headers={"Content-Type": "text/x-fasta"})

        if not req.ok:
            req.raise_for_status()
            logger.error(f"Failed to fetch sequence for {enst}")
        else:
            sequences[enst] = req.text
            with (REFSEQ_CACHE_DIR / f"{enst}.fa").open("w") as f:
                f.write(req.text)

    return sequences


def gb_fetch_from_accession(accession: str) -> SeqRecord:
    """
    Fetch a GenBank record given an accession ID, and cache it locally

    Args:
        accession (str): Accession ID for the record of interest (e.g. NM_023110.3)

    Returns:
        SeqRecord: GenBank record as a SeqRecord object
    """
    cache_file = REFSEQ_CACHE_DIR / f"{accession.replace('.', '_')}.gb"
    if cache_file.exists():
        raise AttributeError("GB Fetch called on a record that is already cached! This should not happen, check your code.")

    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = handle.read()
    with cache_file.open("w") as f:
        f.write(record)
    handle.close()

    return SeqIO.read(cache_file.open("r"), "genbank")


def fetch_exons(accession: str) -> list[SeqRecord]:
    """
    Provide a set of exon sequences for a given accession ID.

    Args:
        accession (str): Accession ID for the transcript of interest (e.g. NM_023110.3)

    Returns:
        list[SeqRecord]: List of exon sequences for the given transcript
    """
    cache_file = REFSEQ_CACHE_DIR / f"{accession.replace('.', '_')}.gb"
    if cache_file.exists():
        with cache_file.open("r") as f:
            record: SeqRecord = SeqIO.read(f, "genbank")
        return extract_exons(record)

    record = gb_fetch_from_accession(accession)
    return extract_exons(record)
