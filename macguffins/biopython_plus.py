"""
    Program | RPINerd, 03/25/26

    Description for program
"""

import logging

from Bio.SeqRecord import Seq, SeqRecord

logging.basicConfig(
    level=logging.DEBUG,
    format=("%(asctime)s | %(levelname)-7s | %(module)-10s | %(lineno)-4d | %(message)s"),
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


def extract_exons(record: SeqRecord) -> list[SeqRecord]:
    """
    Extract exon sequences from a GenBank record

    Args:
        record (SeqRecord): Biopython SeqRecord object containing GenBank data

    Returns:
        list[SeqRecord]: List of SeqRecord objects, each representing an exon sequence
    """
    exons = []
    for feature in record.features:
        if feature.type == "exon":
            exon_seq = feature.extract(record.seq)
            if not isinstance(exon_seq, Seq):
                logger.error(f"Failed to extract exon sequence for {record.id}, {feature}")
                continue
            exon_record = SeqRecord(
                exon_seq,
                id=f"{record.id}_{feature.qualifiers.get('number', ['unknown'])[0]}",
                description=f"Exon {feature.qualifiers.get('number', ['unknown'])[0]} from {record.id}"
            )
            exons.append(exon_record)
    return exons
