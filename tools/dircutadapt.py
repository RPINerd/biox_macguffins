"""
CutAdapt on Directory | RPINerd, 01/18/25

Tiny wrapper to run a cutadapt command, only for a very specific instance at this point
"""

import argparse
import re
import subprocess
from pathlib import Path

from macguffins.macguffin_utils import collect_fastqs


def arg_parse() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Run cutadapt on all fastq files in a directory"
    )

    parser.add_argument(
        "-d",
        "--directory",
        type=Path,
        required=True,
        help="Directory containing fastq files",
    )

    parser.add_argument(
        "-a",
        "--adapters",
        type=Path,
        required=True,
        help="Fasta file of adapter sequences to trim",
    )

    args = parser.parse_args()

    directory = Path(args.directory)
    if not directory.is_dir():
        raise NotADirectoryError(f"{directory} is not a valid directory")
    if not directory.exists():
        raise FileNotFoundError(f"{directory} does not exist")

    return args


def main(args: argparse.Namespace) -> None:
    """Main function to run cutadapt on all fastq files in a directory"""
    directory = Path(args.directory)
    adapters = Path(args.adapters)
    fastqs = collect_fastqs(directory)

    reports_file = Path("reports.tsv")
    if reports_file.exists():
        reports_file.unlink()

    for file in fastqs:
        # TODO set to be universal
        file_info = re.match(
            r"(MMV2-R[0-9]{1,2}_S[0-9]{1,2}_L00[1234]_R)([12])(_001.fastq).gz",
            file.name,
        )

        if not file_info:
            continue

        sample: str = file_info.group(1)
        read: str = file_info.group(2)
        suffix: str = file_info.group(3)

        if read == "2":
            continue

        print(file)
        print(sample + " | " + read + " | " + suffix + "\n")

        paired_read = sample + "2" + suffix + ".gz"
        short_r1 = "short_" + sample + str(read) + suffix
        short_r2 = "short_" + sample + "2" + suffix
        trimmed_r1 = "trimmed_" + sample + str(read) + suffix
        trimmed_r2 = "trimmed_" + sample + "2" + suffix

        cut_cmd = [
            "cutadapt",
            "-j",
            "12",
            "--report=minimal",
            "-n",
            "3",
            "-m",
            "25",
            "-a",
            f"file:{adapters}",
            "-A",
            f"file:{adapters}",
            "--too-short-output",
            short_r1,
            "--too-short-paired-output",
            short_r2,
            "--pair-filter=first",
            "-o",
            trimmed_r1,
            "-p",
            trimmed_r2,
            str(file),
            paired_read,
        ]

        print("cutadapt:\n" + " ".join(cut_cmd) + "\n")

        with reports_file.open("a", encoding="utf-8") as reports:
            reports.write(sample + "\n")
            subprocess.run(cut_cmd, stdout=reports, check=True)


if __name__ == "__main__":
    args = arg_parse()

    main(args)
