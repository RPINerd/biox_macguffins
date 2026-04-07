"""
    Simple Utilities | RPINerd, 03/20/24

    Collection of very basic utilities for use in other scripts.
"""
import gzip
import logging
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from Bio import SeqIO

if TYPE_CHECKING:
    from collections.abc import Generator, Iterator

    from Bio.SeqRecord import SeqRecord

RNA_TRANSLATE = str.maketrans("AUGCaugc", "TACGtacg")
RNA_CONVERT = str.maketrans("Uu", "Tt")
DNA_TRANSLATE = str.maketrans("ATCGatcg", "TAGCtagc")
DNAN_TRANSLATE = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
IUPAC_TRANSLATE = str.maketrans(
    "ATGCRYMKSWHBVDNatgcrymkswhbvdn",
    "TACGYRKMSWDVBHNtacgyrkmswdvbhn"
)
WOBBLE_BASES = {
    "R": ["A", "G"],            # puRine
    "Y": ["C", "T"],            # pYrimidine
    "M": ["A", "C"],            # aMino
    "K": ["G", "T"],            # Keto
    "S": ["G", "C"],            # Strong
    "W": ["A", "T"],            # Weak
    "H": ["A", "C", "T"],       # not G
    "B": ["G", "T", "C"],       # not A
    "V": ["G", "A", "C"],       # not T
    "D": ["G", "A", "T"],       # not C
    "N": ["G", "A", "C", "T"],  # aNy
}

logging.basicConfig(
    level=logging.DEBUG,
    format=("%(asctime)s | %(levelname)-7s | %(module)-10s | %(lineno)-4d | %(message)s"),
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


def collect_fastqs(directory: Path) -> list[Path]:
    """
    Collect all fastq and fastq.gz files in a directory and return list of file paths

    Args:
        directory (Path): Directory to search

    Returns:
        list[Path]: List of fastq and fastq.gz files
    """
    fastq_files = []
    for file in directory.iterdir():
        if file.suffix in {".fastq", ".fq"} or file.suffixes == [".fastq", ".gz"] or file.suffixes == [".fq", ".gz"]:
            fastq_files.append(file)

    if not fastq_files:
        raise FileNotFoundError(f"No fastq files found in {directory}")

    return fastq_files


def collect_fastq_pairs(directory: Path) -> dict[str, tuple[Path, Path]]:
    """
    Collect all fastq and fastq.gz files in a directory and return as a dictionary of samples with paired-end files

    Args:
        directory (Path): Directory to search

    Returns:
        dict[str, tuple[Path, Path]]: Dictionary mapping sample names to tuples of paired-end fastq files

    Raises:
        ValueError: If any pairs are missing a read
    """
    fastq_files = collect_fastqs(directory)
    sample_dict: dict[str, list[Path | None]] = {}
    for file in fastq_files:
        sample_name, read_num, _ = extract_sample_info(file)
        if sample_name not in sample_dict:
            sample_dict[sample_name] = [None, None]
        if read_num == 1:
            sample_dict[sample_name][0] = file
        elif read_num == 2:
            sample_dict[sample_name][1] = file

    # Convert lists to tuples and check for missing pairs
    for sample_name, files in sample_dict.items():
        if None in files:
            raise ValueError(f"Missing pair for sample {sample_name}: {files}")
        sample_dict[sample_name] = tuple(files)

    return sample_dict


def contains_n_consecutive(n: int, lst: list, sort: bool = False) -> bool:
    """
    Check if an integer list contains n or more consecutive numbers

    e.g. n = 3, lst = [1, 2, 3, 6, 10]
        returns True because list contains 3 consecutive numbers (1,2,3)
        n = 4, lst = [1, 4, 5, 6, 10]
        returns False because the longest sequence of consescutive numbers is only 3  - (4,5,6)

    Args:
        n (int): Number of consecutive numbers to check for
        lst (list): List of integers to check
        sort (bool): Sort the list before checking

    Returns:
        bool: True if the list contains n or more consecutive numbers
    """
    if sort:
        lst = sorted(lst)

    prev = lst[0]
    count = 1
    for idx, e in enumerate(lst):
        if e - prev == 1:
            count += 1
        else:
            count = 1
        if count == n:
            return True
        prev = e

    return False


def extract_sample_info(filename: Path) -> tuple[str, int, int]:
    """
    Parse a file name to derive info about the samplename, read pair and lane

    Function expects file name to be in reasonable format with different designators separated by underscores.

    ABCRUN123_L001_R1_001.fastq.gz
    DE-F-456_R2_L004.fastq

    Args:
        filename (Path): Input file path

    Returns:
        Sample Info (tuple[str, int, int])
        - Sample ID
        - Read pair number
        - Lane number
    Raises:
        ValueError: If the provided filename has spaces in it
    """
    if " " in filename.name:
        raise ValueError(f"I will not work with files that contain spaces, on principal. ({filename=})")

    sample_name = ""
    read_num = 0
    lane_num = 1

    name_witout_suffixes = filename.stem
    for suffix in filename.suffixes:
        name_witout_suffixes = name_witout_suffixes.removesuffix(suffix)
    components = name_witout_suffixes.split("_")

    logger.debug(f"Extracting sample info from {filename.name} with components {components}")
    sample_name = components[0]
    for comp in components[1:]:
        if comp.startswith(("l", "L")):
            lane_num = int(comp.strip("lL"))
        elif comp.startswith(("R", "r")):
            read_num = int(comp.strip("Rr"))

    if lane_num not in range(1, 9):
        raise ValueError(f"An unexpected lane number was detected for {filename}: {lane_num}")
    if read_num not in {1, 2}:
        raise ValueError(f"Weird read number ({read_num}) for {filename}. Reads must be either 1 or 2")

    return sample_name, read_num, lane_num


def look_backward_match(iterable: list | tuple, start: int, char: str) -> int:
    """
    Look behind in an iterable for the next point where a character is the same

    Args:
        iterable (list | tuple): A list/tuple to look backward through
        start (int): The initial index to being from
        char (str): The character to look for the end of in the sequence

    Returns:
        int: The index of the next point where the character is the same

    Raises:
        ValueError: If no match is found
    """
    idx = start
    end_idx = None
    while not end_idx and idx > 0:
        idx -= 1
        if iterable[idx] == char:
            return idx + 1

    raise ValueError(f"No match found looking backwards from index {start} along interable:\n{iterable}")


def look_backward_miss(iterable: list | tuple, start: int, char: str) -> int:
    """
    Look behind in an iterable for the next point where a character is different

    Args:
        iterable (list | tuple): A list/tuple to look backward through
        start (int): The initial index to being from
        char (str): The character to look for the end of in the sequence

    Returns:
        int: The index of the next point where the character is different

    Raises:
        ValueError: If no match is found
    """
    idx = start
    end_idx = None
    while not end_idx and idx > 0:
        idx -= 1
        if iterable[idx] != char:
            end_idx = idx

    if end_idx is None:
        raise ValueError(
            f"No match found looking backwards from index {start} along interable:\n{iterable[start:len(iterable)]}"
        )
    return end_idx


def look_forward_match(iterable: list | tuple, start: int, char: str) -> int:
    """
    Look ahead in an iterable for the next point where a character is the same

    Args:
        iterable (list | tuple): A list/tuple to look forward through
        start (int): The initial index to being from
        char (str): The character to look for the end of in the sequence

    Returns:
        int: The index of the next point where the character is the same

    Raises:
        ValueError: If no match is found
    """
    idx = start
    end_idx = None
    while not end_idx and idx < len(iterable):
        idx += 1
        if iterable[idx] == char:
            return idx - 1

    raise ValueError(
        f"No match found looking forwards from index {start} along interable:\n{iterable[start:len(iterable)]}"
    )


def look_forward_miss(iterable: list | tuple, start: int, char: str) -> int:
    """
    Look ahead in an iterable for the next point where a character is different

    Args:
        iterable (list | tuple): A list/tuple to look forward through
        start (int): The initial index to being from
        char (str): The character to look for the end of in the sequence

    Returns:
        int: The index of the next point where the character is different

    Raises:
        ValueError: If no match is found
    """
    idx = start
    end_idx = None
    while not end_idx and idx < len(iterable):
        idx += 1
        if iterable[idx] != char:
            end_idx = idx

    if end_idx is None:
        raise ValueError(
            f"No match found looking forwards from index {start} along interable:\n{iterable[start:len(iterable)]}"
        )
    return end_idx


def ret_idt_repr(seq: str) -> str:
    """
    Return IDT representation of an oligodesign2 sequence

    Uppercases all bases and adds a + before LNA bases (upper cased bases in OD2)

    Args:
        seq (str): Oligodesign2 sequence with/without LNAs

    Returns:
        str: IDT sequence
    """
    # Already in IDT format
    if seq.find("+") != -1:
        print("Sequence already in IDT format")
        return seq
    idt_seq = []
    for alphabet in seq:
        assert alphabet.lower() in {"a", "c", "g", "t"}, seq
        # Uppercase Base - LNA
        if alphabet.upper() == alphabet:
            alphabet = "+" + alphabet
        # Lowercase Base
        else:
            alphabet = alphabet.upper()
        idt_seq.append(alphabet)
    return "".join(idt_seq)


def ret_od2_repr(seq: str) -> str:
    """
    Return Oligodesign2 representation of an IDT sequence

    Args:
        seq (str): IDT sequence with/without LNAs

    Returns:
        str: OD2 sequence
    """
    is_lna = False
    od2_seq = []
    for i, a in enumerate(seq):
        assert a.upper() == a, "IDT bases should be upper case"
        # Next base is an LNA base
        if a == "+":
            is_lna = True
            continue
        if is_lna:
            od2_seq.append(a)
        else:
            od2_seq.append(a.lower())
        is_lna = False

    return "".join(od2_seq)


def revcomp(seq: str) -> str:
    """
    Reverse complement a sequence

    Args:
        seq (str): Genetic sequence

    Returns:
        str: Reverse complemented sequence
    """
    return seq.translate(DNAN_TRANSLATE)[::-1]


def smart_iter(filepath: Path) -> Iterator[str]:
    """
    Iterate over a file, handling gzip files automatically.

    Args:
        filepath (Path): Path to the file.

    Returns:
        Iterator over lines in the file.
    """
    if filepath.suffix == ".gz":
        with gzip.open(filepath, "r") as f:
            for line in f:
                yield line.decode("utf-8").strip("\n")
    else:
        with Path.open(filepath) as f:
            for line in f:
                yield line.strip("\n")


def smart_open(filepath: Path) -> TextIO:
    """
    Open a file, handling gzip files automatically.

    Args:
        filepath (Path): Path to the file.

    Returns:
        File object.
    """
    if filepath.suffix == ".gz":
        return gzip.open(filepath, "rt")
    return filepath.open("r", encoding="utf-8")


def smart_fq_zip(file_a: Path, file_b: Path) -> Generator[tuple[SeqRecord, SeqRecord]]:
    """
    Zip-like function that yields iterators of SeqRecords from two input fastq files

    Args:
        file_a (Path): Fastq file 1
        file_b (Path): Fastq file 2

    Yields:
        Tuple of 2 Biopython SeqRecords
    """
    records_a: Generator[SeqRecord] = SeqIO.parse(smart_open(file_a), "fastq")
    records_b: Generator[SeqRecord] = SeqIO.parse(smart_open(file_b), "fastq")
    yield from zip(records_a, records_b)
