"""
Extract index sequences from paired-end FASTQ files (R1/R2) and write to I1/I2 FASTQ files.

Usage:
    python fq_index_extract.py -r R1.fastq.gz
"""

import argparse
import gzip
import logging
from pathlib import Path
from typing import IO, Any

from Bio import SeqIO


def open_fastq(path: Path, mode: str = "rt") -> IO[Any] | gzip.GzipFile:
    """Open a FASTQ file, handling gzip if needed."""
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return Path.open(path, mode)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Extract index sequences from paired FASTQ files.")
    parser.add_argument("-r", "--read_file", type=Path, required=True, help="Path to R1 FASTQ file (can be .gz)")
    parser.add_argument("-1", "--index_1", type=Path, required=False, help="Output path for I1 FASTQ file")
    parser.add_argument("-2", "--index_2", type=Path, required=False, help="Output path for I2 FASTQ file")
    return parser.parse_args()


def main(read_file: Path, i1_file: Path, i2_file: Path) -> None:
    """
    Extract index sequences from paired R1/R2 FASTQ files and write to I1/I2 files.

    Args:
        r1_path: Path to R1 FASTQ file
        i1_path: Output path for I1 FASTQ file
        i2_path: Output path for I2 FASTQ file
    """
    with open_fastq(read_file) as r1_handle, \
        open_fastq(i1_file, "wt") as i1_handle, open_fastq(i2_file, "wt") as i2_handle:
        for rec in SeqIO.parse(r1_handle, "fastq"):
            # Extract index from read id (assumes Illumina format: ...:INDEX1+INDEX2)
            rec_id_parts: list[str] = rec.id.split(":")
            if len(rec_id_parts) < 2:
                logging.warning(f"Malformed read ID: {rec.id}")
                continue
            i1_seq = rec_id_parts[-1].split("+")[0]
            i2_seq = rec_id_parts[-1].split("+")[1]
            # Write to I1 and I2 as FASTQ (dummy quality)
            i1_handle.write(f"@{rec.id}\n{i1_seq}\n+\n{'I' * len(i1_seq)}\n")
            i2_handle.write(f"@{rec.id}\n{i2_seq}\n+\n{'I' * len(i2_seq)}\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()
    read_file = Path(args.read_file)
    if not args.index_1 or not args.index_2:
        read_file_name = read_file.name
        i_file_name_1 = read_file_name.replace("_R1", "_I1")
        i_file_name_2 = read_file_name.replace("_R1", "_I2")
        index_file_1 = read_file.parent / i_file_name_1
        index_file_2 = read_file.parent / i_file_name_2
    else:
        index_file_1 = Path(args.index_1)
        index_file_2 = Path(args.index_2)
    logging.info(f"Extracting indices from {read_file} to {index_file_1} and {index_file_2}")
    assert index_file_1 != index_file_2, "I1 and I2 output files must be different!"
    assert not index_file_1.exists() and not index_file_2.exists(), "I1/I2 output files must not already exist!"
    main(read_file, index_file_1, index_file_2)
    logging.info(f"Index extraction complete: {index_file_1}, {index_file_2}")
