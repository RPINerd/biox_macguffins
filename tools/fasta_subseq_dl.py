"""
    Fasta Subsequence Download Tool | RPINerd, 11/26/24

    Given a file with Accession IDs and start/stop coordinates, download the fasta sequence for each entry.
    Entry format is as follows:
        CP060352.1:280230-280675
"""

import argparse
import logging
from io import TextIOWrapper
from pathlib import Path
from time import sleep

from Bio import Entrez

logging.basicConfig(
    level=logging.DEBUG,
    format=("%(asctime)s | %(levelname)-7s | %(module)-10s | %(lineno)-4d | %(message)s"),
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, help="Input file with Accession IDs and start/stop coordinates")
    parser.add_argument("-b", "--bed", action="store_true", help="Flag for bed format input", required=False, default=False)
    parser.add_argument("-o", "--output", type=Path, help="Directory to output file(s) to.", required=False)
    args = parser.parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input file {args.input} was not found!")

    if not args.output:
        args.output = args.input.parent
        logger.debug(f"No output dir specified, using input parent: {args.output.as_posix()}")

    return args


def coord_extract(lines: TextIOWrapper) -> dict[str, list[tuple[str, str]]]:
    """Extract ID and coordinates from standard format"""
    logger.debug("Extracting coordinates from traditional ID list...")

    regions: dict[str, list[tuple[str, str]]] = {}
    total_regions: int = 0
    for line in lines:
        enst, coords = line.strip().split(":")
        start, stop = coords.split("-")
        try:
            regions[enst].append((start, stop))
        except KeyError:
            regions[enst] = [(start, stop)]
        except Exception as e:
            raise Exception(f"Something went wrong during coordinate extraction:\n{e}")
        total_regions += 1

    logger.debug(f"Extracted {total_regions=} from file...")
    return regions


def bed_extract(lines: TextIOWrapper) -> dict[str, list[tuple[str, str]]]:
    """Extract ID and coordinates from bed format"""
    logger.debug("Extracting coordinates from BED-format ID list...")

    regions: dict[str, list[tuple[str, str]]] = {}
    total_regions: int = 0
    for line in lines:
        enst, start, stop = line.strip().split()
        try:
            regions[enst].append((start, stop))
        except KeyError:
            regions[enst] = [(start, stop)]
        except Exception as e:
            raise Exception(f"Something went wrong during coordinate extraction:\n{e}")
        total_regions += 1

    logger.debug(f"Extracted {total_regions=} from file...")
    return regions


def main(input_file: Path, bed_format: bool, output_location: Path) -> None:
    """Main function"""
    logger.info(f"Pulling fasta sequences from regions file {input_file.name}")

    with input_file.open("r", encoding="utf-8") as id_list:
        regions_to_download = {}
        if bed_format:
            regions_to_download = bed_extract(id_list)
        else:
            regions_to_download = coord_extract(id_list)
    logger.info(f"Extracted {len(regions_to_download)} accessions...")

    for enst, regions in regions_to_download.items():
        logger.debug(f"Fetching {enst}...")
        handle = Entrez.efetch(db="nucleotide", id=enst, rettype="fasta")
        record = handle.read()
        handle.close()
        lines = record.split("\n")
        header = lines[0]
        sequence = "".join(lines[1:])

        outfile = output_location / f"{enst.split(".")[0]}.fasta"
        logger.debug(f"Writing {len(regions)} regions to {outfile.as_posix()}...")
        with outfile.open("w", encoding="utf-8") as out:
            for start, stop in regions:
                out.write(f"{header}\n{sequence[int(start):int(stop)]}\n")
        sleep(3)


if __name__ == "__main__":
    args = parse_args()
    main(args.input, args.bed, args.output)
