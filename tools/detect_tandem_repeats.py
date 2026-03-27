"""
    Tandem Repeat Inspector | RPINerd, 03/26/26

    Given an assembled FASTQ read set (R1/R2 -> contigs), identify reads that contain tandem repeats of a specified gene.

    Usage:
        python detect_tandem_repeats.py -r contigs.fasta -a gene_accession -o output.txt
"""

import argparse
import logging
import sys
from pathlib import Path

import edlib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
from macguffins.macguffin_acquire import fetch_exons

logging.basicConfig(
    level=logging.DEBUG,
    format=("%(asctime)s | %(levelname)-7s | %(module)-10s | %(lineno)-4d | %(message)s"),
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)

MIN_EXON_SIMILARITY = 0.90


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--read_file",
        type=Path,
        required=True,
        help="Path to assembled reads in FASTA format",
    )
    parser.add_argument(
        "-a",
        "--accession",
        type=str,
        required=True,
        help="Accession ID of the gene to search for tandem repeats",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Path to output file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    return args


def find_exon_hits(contig: SeqRecord, exons: list[SeqRecord]) -> list[str]:
    """
    Check for exon hits within the contig sequence and log positions and exon numbers if found.

    Args:
        contig: A SeqRecord object representing the contig sequence.
        exons: A list of exons to search for.

    Returns:
        An array of exon numbers that were found, in the order they appear in the contig.
    """
    logger.debug(f"Scanning contig {contig.id} for exon hits...")
    contig_seq = str(contig.seq)
    hits = []
    for exon in exons:
        exon_id = exon.id or ""
        exon_num = exon_id.split("_")[-1] if exon_id else "unknown"
        exon_seq = str(exon.seq)
        max_edits = int((1.0 - MIN_EXON_SIMILARITY) * len(exon_seq))
        result = edlib.align(exon_seq, contig_seq, mode="HW", task="locations", k=max_edits)
        if result["editDistance"] == -1:
            continue

        similarity = 1.0 - (result["editDistance"] / len(exon_seq))
        if similarity >= MIN_EXON_SIMILARITY:
            start_pos = result["locations"][0][0]
            logger.debug(
                "Found exon %s in contig %s at position %s with similarity %.3f",
                exon_num,
                contig.id,
                start_pos,
                similarity,
            )
            hits.append(exon_num)
    return hits


def main(args: argparse.Namespace) -> None:
    """Main function to detect tandem repeats in assembled reads"""
    logger.info(f"Starting tandem repeat detection for {args.accession} in {args.read_file}")
    exons = fetch_exons(args.accession)
    logger.info(f"Fetched {len(exons)} exons for {args.accession}")

    if logger.isEnabledFor(logging.DEBUG):
        for exon in exons:
            logger.debug(f"Exon {exon.id}: {len(exon.seq)} bp")

    contigs_generator = SeqIO.parse(args.read_file, "fasta")
    exon_maps = {}
    logger.info("Scanning contigs for exon hits...")
    for contig in contigs_generator:
        hits = find_exon_hits(contig, exons)
        if hits:
            exon_maps[contig.id] = hits

    if args.output:
        with args.output.open("w") as f:
            for contig_id, hits in exon_maps.items():
                f.write(f"{contig_id}\t{','.join(hits)}\n")
        logger.info(f"Wrote results to {args.output}")
    else:
        for contig_id, hits in exon_maps.items():
            print(f"{contig_id}\t{','.join(hits)}")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
