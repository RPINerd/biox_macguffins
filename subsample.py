"""
    SeqTK Subsampling | RPINerd, 05/21/25

    Wrapper around SeqTK's sample function to generate a subset fastq file based
    on the file given by the user and the number of subsamples requested.
"""

import argparse
import logging
import os
import random
import re
from pathlib import Path

logger = logging.getLogger(__name__)


# TODO allow parsing of a file to bulk subsample
def parse_tsv(file: Path) -> list[Path]:
    """
    Parse a TSV file to extract the paths of fastq files

    Args:
        file (Path): Path to the input TSV file

    Returns:
        list[Path]: List of paths to fastq files
    """
    workdir = file.parent
    with Path.open(file, "r") as f:
        lines = f.readlines()
        input_files = []
        for i, line in enumerate(lines):
            rawline = line.strip()

            if rawline.startswith("#") or not rawline:
                continue

            cols = rawline.split("\t")
            if len(cols) > 1:
                logger.debug(f"Multiple columns found on line {i}: {rawline}.. Ignoring extra columns")

            input_file = workdir.joinpath(cols[0])

            if not input_file.exists():
                logger.warning(f"Input file {input_file} does not exist! Skipping...")
                continue

            input_files.append(input_file)

    logger.debug(f"Parsed {len(input_files)} input files from TSV")
    return input_files


def arg_parse() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate a fastq file with a random subset of reads from the given input file"
    )
    parser.add_argument("-f", "--file", type=Path, required=True, help="File, or tsv file containing reads you wish to subset")
    parser.add_argument("-v", "--verbose", type=bool, action="store_true", help="Lots of status messages")
    parser.add_argument(
        "-n", "--number", type=int, default=2000, help="Number of reads desired in the subset (default: 2000)"
    )
    return parser.parse_args()


def main(file: Path, sub_size: int) -> None:
    """
    Main function to handle the subsampling of fastq files

    Args:
        file (Path): Path to the input file
        outfile (str): Path to the output file
        sub_size (int): Number of reads to subsample
    """
    if not file.exists():
        logger.error("Input file does not exist!")
        raise FileNotFoundError(f"Input file {file} does not exist!")

    input_files: list[Path] = []
    if file.suffix == ".tsv":
        logger.debug("Input file is a tsv file")
        # TODO implement TSV parsing
        input_files = parse_tsv(file)
    elif file.suffix == ".fastq.gz":
        logger.debug("Input file is a fastq.gz file")
        input_files = [file]
    elif file.suffix == ".fastq":
        logger.debug("Input file is a fastq file")
        input_files = [file]
    else:
        logger.error("Input file is not a valid fastq or tsv file")
        raise ValueError(f"Input file {file} is not a valid fastq or tsv file")

    seed = random.randint(1, 999)
    logger.debug(f"Seed value for this run is: {seed}")

    for input_file in input_files:
        logger.debug(f"Processing input file: {input_file}")
        file_reg = re.search(r"(^.*)_R[12](.*fastq.gz)", input_file.name)
        if not file_reg:
            logger.error("Input file name does not match expected pattern")
            raise ValueError(f"Input file {input_file} does not match expected pattern")
        sample_id = file_reg.group(1)
        suffix = file_reg.group(2)
        read1_file = f"{sample_id}_R1{suffix}"
        read2_file = f"{sample_id}_R2{suffix}"
        logger.debug(f"Sample ID: {sample_id}")
        logger.debug(f"File Suffix: {suffix}")
        logger.debug(f"Read files: {read1_file}, {read2_file}")

        # Prepare output files
        outfile_r1 = f"{sample_id}_R1.{sub_size}.fastq"
        outfile_r2 = f"{sample_id}_R2.{sub_size}.fastq"
        logger.debug(f"Writing out to files: {outfile_r1}, {outfile_r2}")

        # Execute the calls to seqtk
        sub1 = f"seqtk sample -s {seed} {read1_file} {sub_size} > {outfile_r1}"
        logger.debug("Read 1 call generated and executing: \n" + sub1)
        os.system(sub1)
        sub2 = f"seqtk sample -s {seed} {read2_file} {sub_size} > {outfile_r2}"
        logger.debug("Read 2 call generated and executing: \n" + sub2)
        os.system(sub2)


if __name__ == "__main__":

    # Verify that seqtk is installed
    if os.system("seqtk --version") != 0:
        logger.error("SeqTK is not installed. Please install it to use this script.")
        raise OSError("SeqTK is not installed. Please install it to use this script.")

    args = arg_parse()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    main(args.file, args.number)
