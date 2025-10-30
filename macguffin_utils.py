"""
    Simple Utilities | RPINerd, 03/20/24

    Collection of very basic utilities for use in other scripts.
"""
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
