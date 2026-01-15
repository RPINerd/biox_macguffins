"""
    Simple Utilities | RPINerd, 03/20/24

    Collection of very basic utilities for use in other scripts.
"""
import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

RNA_TRANSLATE = str.maketrans("AUGCaugc", "TACGtacg")
RNA_CONVERT = str.maketrans("Uu", "Tt")
DNA_TRANSLATE = str.maketrans("ATCGatcg", "TAGCtagc")
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


def collect_fastqs(directory: Path) -> list[Path]:
    """
    Collect all fastq and fastq.gz files in a directory and return as a dictionary of samples

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


def extract_sample_info(filename: str) -> tuple[str, int, int]:
    """
    Parse a file name to derive info about the samplename, read pair and lane

    Function expects file name to be in reasonable format with different designators separated by underscores.

    ABCRUN123_L001_R1_001.fastq.gz
    DE-F-456_R2_L004.fastq

    Args:
        filename (str): Name of input file

    Returns:
        Sample Info (tuple[str, int, int])
        - Sample ID
        - Read pair number
        - Lane number
    Raises:
        ValueError: If the provided filename has spaces in it
    """
    if " " in filename:
        raise ValueError(f"I will not work with files that contain spaces, on principal. ({filename=})")

    sample_name = ""
    read_num = 0
    lane_num = 0

    components = filename.split("_")
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
