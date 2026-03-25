"""
Brute force script to detect fusion reads in paired-end FASTQ files.

Fusions are extracted from 2 sets of FASTA files, "fusionA" and "fusionB"

1. Build a k-mer index from each reference FASTA(s).
2. Remove k-mers that are shared between the two genes to eliminate
    false-positive hits from homologous sequence.
3. Scan every read in R1 and R2 for evidence of:
    Spanning reads: a single read contains k-mers from both fusA and fusB - the read likely crosses the fusion junction.
    Chimeric pairs: one mate maps to fusA and the other to fusB - the fragment bridges the junction but neither read individually spans it.

Output:
    Prints total raw reads, spanning-read count, and chimeric-pair count.

Usage:
    uv run detect_fusions.py --fusionA first_fusion_partner.fasta --fusionB second_fusion_partner.fasta --fastq-dir fastqs/
"""

import argparse
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from macguffins.macguffin_utils import collect_fastq_pairs, revcomp, smart_fq_zip

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

DEFAULT_KMER_SIZE = 20
_DNA_BASES = ("A", "C", "G", "T")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fusionA",
        "-a",
        required=True,
        nargs="+",
        type=Path,
        metavar="FASTA",
        help="One or more reference FASTA files of your first fusion sequence.",
    )
    parser.add_argument(
        "--fusionB",
        "-b",
        required=True,
        nargs="+",
        type=Path,
        metavar="FASTA",
        help="One or more reference FASTA files of your second fusion sequence.",
    )
    parser.add_argument(
        "--fastq-dir",
        "-i",
        required=True,
        type=Path,
        metavar="DIR",
        help="Directory containing *_R1.fastq.gz and *_R2.fastq.gz files.",
    )
    parser.add_argument(
        "--kmer-size",
        "-k",
        type=int,
        default=DEFAULT_KMER_SIZE,
        metavar="K",
        help=f"K-mer length used for matching (default: {DEFAULT_KMER_SIZE}).",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        choices=(0, 1),
        default=0,
        metavar="N",
        help="Maximum mismatches allowed per k-mer match (default: 1).",
    )
    parser.add_argument(
        "--out-dir",
        "-o",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "Directory to write hit reads as gzipped FASTQs. "
            "For each input pair, two files are written: "
            "<stem>_fusion_R1.fastq.gz and <stem>_fusion_R2.fastq.gz. "
            "Omit to skip writing output FASTQs."
        ),
    )
    args = parser.parse_args()

    if args.out_dir is None:
        args.out_dir = args.fastq_dir.parent / "fusion_hits"
    args.out_dir.mkdir(parents=True, exist_ok=True)

    return parser.parse_args()


def build_kmer_set(fastas: list[Path], k: int) -> set[str]:
    """
    Build a set of k-mers (forward + reverse complement) from fasta file(s).

    K-mers containing 'N' are silently skipped.

    Args:
        fastas (list[Path]): One or more referece fastas to build kmers from

    Returns:
        Set of k-mer strings.
    """
    kmers: set[str] = set()
    for fasta in fastas:
        for record in SeqIO.parse(fasta, "fasta"):
            n = len(record.seq)
            for i in range(n - k + 1):
                kmer = record.seq[i : i + k]
                if "N" not in kmer:
                    kmers.add(kmer)
                    kmers.add(revcomp(kmer))
    return kmers


def read_has_kmers(
    sequence: str,
    kmer_set: set[str],
    k: int,
    max_mismatches: int,
) -> bool:
    """
    Return True if {sequence} contains at least one k-mer in {kmer_set}.

    Args:
        sequence: Raw read sequence (any case).
        kmer_set: Pre-built k-mer set (includes reverse complements).
        k: K-mer length.
        max_mismatches: Maximum mismatches allowed between a read window and
            a k-mer from kmer_set.

    Returns:
        True if any k-mer from the read is present in kmer_set.
    """
    seq = sequence.upper()
    n = len(seq)
    for i in range(n - k + 1):
        window = seq[i : i + k]
        if "N" in window:
            continue

        # Fast path: exact match
        if window in kmer_set:
            return True

        if max_mismatches < 1:
            continue

        # One-mismatch path: try all single-base substitutions.
        # TODO replace with more performant fuzzy matching
        for j, original_base in enumerate(window):
            prefix = window[:j]
            suffix = window[j + 1 :]
            for new_base in _DNA_BASES:
                if new_base == original_base:
                    continue
                candidate = prefix + new_base + suffix
                if candidate in kmer_set:
                    return True
    return False


def scan_fastq_pair(
    r1_path: Path,
    r2_path: Path,
    fusion_a_kmers: set[str],
    fusion_b_kmers: set[str],
    k: int,
    max_mismatches: int,
    out_path: Path | None = None,
) -> dict[str, int]:
    """
    Scan a FASTQ pair and count fusion evidence.

    total_reads - every individual read (R1 + R2).
    spanning_reads - individual reads that contain k-mers from both
        fusion sets (most likely junction-spanning reads).
    chimeric_pairs - read pairs where one mate hits FusionA and the other
        hits FusionB, but neither is individually spanning.
    any_fusion_pairs - read pairs with any combination of the above
        (spanning or chimeric).

    Args:
        r1_path: Path to R1 FASTQ.
        r2_path: Path to R2 FASTQ.
        fusion_a_kmers: FusionA-specific k-mer set.
        fusion_b_kmers: FusionB-specific k-mer set.
        k: K-mer length.
        max_mismatches: Maximum mismatches allowed when matching k-mers.
        out_path: Optional output path for hit reads (gzip-compressed).

    Returns:
        Dictionary with counts: total_reads, spanning_reads,
        chimeric_pairs, any_fusion_pairs.
    """
    total_reads = 0
    spanning_reads = 0
    chimeric_pairs = 0
    any_fusion_pairs = 0
    fusion_hits: list[tuple[SeqRecord, SeqRecord]] = []

    for record_1, record_2 in smart_fq_zip(r1_path, r2_path):
        total_reads += 2  # one R1 + one R2 per pair

        r1_fus_1 = read_has_kmers(record_1.seq, fusion_a_kmers, k, max_mismatches)
        r1_fus_2 = read_has_kmers(record_1.seq, fusion_b_kmers, k, max_mismatches)
        r2_fus_1 = read_has_kmers(record_2.seq, fusion_a_kmers, k, max_mismatches)
        r2_fus_2 = read_has_kmers(record_2.seq, fusion_b_kmers, k, max_mismatches)

        # Spanning: individual read holds k-mers from both genes
        r1_spanning = r1_fus_1 and r1_fus_2
        r2_spanning = r2_fus_1 and r2_fus_2

        if r1_spanning:
            spanning_reads += 1
        if r2_spanning:
            spanning_reads += 1

        # Chimeric pair: mates hit different genes (but neither spans alone)
        pair_chimeric = (
            not r1_spanning
            and not r2_spanning
            and ((r1_fus_1 and r2_fus_2) or (r1_fus_2 and r2_fus_1))
        )
        if pair_chimeric:
            chimeric_pairs += 1

        is_hit = r1_spanning or r2_spanning or pair_chimeric
        if is_hit:
            any_fusion_pairs += 1
            fusion_hits.append((record_1, record_2))

    if not out_path:
        out_path = r1_path.parent

    r1_fus_file = (out_path / f"{r1_path.name}.fusion.fastq").open("wt")
    r2_fus_file = (out_path / f"{r2_path.name}.fusion.fastq").open("wt")

    for fus_1, fus_2 in fusion_hits:
        SeqIO.write(fus_1, r1_fus_file, "fastq")
        SeqIO.write(fus_2, r2_fus_file, "fastq")

    r1_fus_file.close()
    r2_fus_file.close()

    return {
        "total_reads": total_reads,
        "spanning_reads": spanning_reads,
        "chimeric_pairs": chimeric_pairs,
        "any_fusion_pairs": any_fusion_pairs,
    }


def write_report(args: argparse.Namespace, totals: dict[str, int]) -> None:
    """Print out stats for the run"""
    total = totals["total_reads"]
    spanning = totals["spanning_reads"]
    chimeric = totals["chimeric_pairs"]
    any_fusion = totals["any_fusion_pairs"]
    fusion_a_names = ", ".join(p.name for p in args.fusionA)
    fusion_b_names = ", ".join(p.name for p in args.fusionB)

    logger.info("=" * 55)
    logger.info("Fusion Report Follows:")
    report = f"""K-mer size              : {args.kmer_size}\n
        Max mismatches          : {args.max_mismatches}\n
        Reference Fusion A      : {fusion_a_names}\n
        Reference Fusion B      : {fusion_b_names}\n
        Output directory        : {args.out_dir}\n{"-" * 55}\n
        Total raw reads         : {total:>12,}\n
        Total read pairs        : {total // 2:>12,}\n{"-" * 55}\n
        Spanning reads          : {spanning:>12,} ({spanning / total * 100:.4f}% of reads)\n
        Chimeric pairs          : {chimeric:>12,} ({chimeric / (total // 2) * 100:.4f}% of pairs)\n
        Pairs with any evidence : {any_fusion:>12,} ({any_fusion / (total // 2) * 100:.4f}% of pairs)"\n\n
        Definitions:\n
            \tSpanning reads  — individual reads containing k-mers from\n
            \t\tBOTH fusions (likely junction-spanning).\n
            \tChimeric pairs  — read pairs where R1 hits one gene and R2\n
            \t\thits the other (fragment spans the junction).\n
            \tPairs w/ any    — union of spanning and chimeric evidence.\n
    """
    logger.info("=" * 55 + "\n" + report)


def main(args: argparse.Namespace) -> None:
    """Driver for tool. Build k-mer indices, scan FASTQs, report."""
    logger.info(f"Parsing FusionA reference(s) into k-mers: {args.fusionA}")
    fusion_a_kmers = build_kmer_set(args.fusionA, args.kmer_size)
    logger.info(f"  FusionA unique k-mers (before dedup): {len(fusion_a_kmers)}")

    logger.info(f"Parsing FusionB reference(s) into k-mers: {args.fusionB}")
    fusion_b_kmers = build_kmer_set(args.fusionB, args.kmer_size)
    logger.info(f"  FusionB unique k-mers (before dedup): {len(fusion_b_kmers)}")

    shared = fusion_a_kmers & fusion_b_kmers
    if shared:
        logger.info(
            f"Removing {len(shared)} k-mers shared between fusions "
            "(would cause false positives)."
        )
        fusion_a_kmers -= shared
        fusion_b_kmers -= shared
    logger.info(
        f"  Final k-mer counts — FusionA: {len(fusion_a_kmers)}  |  FusionB: {len(fusion_b_kmers)}"
    )

    logger.info("Collecting fastq read pair files...")
    fastq_pairs = collect_fastq_pairs(args.fastq_dir)

    totals: dict[str, int] = {
        "total_reads": 0,
        "spanning_reads": 0,
        "chimeric_pairs": 0,
        "any_fusion_pairs": 0,
    }

    for sample, (r1_path, r2_path) in fastq_pairs.items():
        logger.info(f"Processing sample {sample}")
        logger.info(f"    R1: {r1_path.name}")
        logger.info(f"    R2: {r2_path.name}")

        counts = scan_fastq_pair(
            r1_path,
            r2_path,
            fusion_a_kmers,
            fusion_b_kmers,
            args.kmer_size,
            args.max_mismatches,
            out_path=args.out_dir
        )
        for key in totals:
            totals[key] += counts[key]

        n_pairs = counts["total_reads"] // 2
        logger.info(
            f"  Pairs: {n_pairs} | Spanning reads: {counts["spanning_reads"]} | Chimeric pairs: {counts["chimeric_pairs"]}"
        )

    write_report(args, totals)
    logger.info(f"Fastqs with fusion hits are located in {args.out_dir}")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
