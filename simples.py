"""
    Simple Utilities | RPINerd, 03/20/24

    Collection of very basic utilities for use in other scripts.
"""


def contains_n_consecutive(n, lst, sort=False):
    """
    Check if an integer list contains n or more consecutive numbers

    e.g. n = 3, lst = [1, 2, 3, 6, 10]
        returns True because list contains 3 consecutive numbers (1,2,3)
        n = 4, lst = [1, 4, 5, 6, 10]
        returns False because the longest sequence of consescutive numbers is only 3  - (4,5,6)

    :param int: n: Number of consecutive numbers to check for
    :param list: lst: List of integers to check
    :param bool: sort: Sort the list before checking
    :rtype: bool
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


def ret_idt_repr(seq):
    """
    Return IDT representation of an oligodesign2 sequence.
    Uppercases all bases and adds a + before LNA bases (upper cased bases in OD2)

    :param str: seq: Oligodesign2 sequence with/without LNAs
    :rtype: str
    """

    # Already in IDT format
    if seq.find("+") != -1:
        print("Sequence already in IDT format")
        return seq
    idt_seq = []
    for alphabet in seq:
        assert alphabet.lower() in ["a", "c", "g", "t"], seq
        # Uppercase Base - LNA
        if alphabet.upper() == alphabet:
            alphabet = "+" + alphabet
        # Lowercase Base
        else:
            alphabet = alphabet.upper()
        idt_seq.append(alphabet)
    return "".join(idt_seq)


def ret_od2_repr(seq):
    """
    Return Oligodesign2 representation of an IDT sequence

    :param str: seq: IDT sequence
    :rtype: str
    """

    is_LNA = False
    od2_seq = []
    for i, a in enumerate(seq):
        assert a.upper() == a, "IDT bases should be upper case"
        # Next base is an LNA base
        if a == "+":
            is_LNA = True
            continue
        if is_LNA:
            od2_seq.append(a)
        else:
            od2_seq.append(a.lower())
        is_LNA = False

    return "".join(od2_seq)


def revcomp(seq) -> str:
    """
    Return the reverse complement of a sequence

    :param str: seq: Sequence to reverse complement
    :rtype: str
    """

    return seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]


def translateRNA(seq) -> str:
    """
    Return the translation of a RNA sequence

    :param str: seq: RNA sequence to translate
    :rtype: str
    """

    return seq.translate(str.maketrans("AUGCaugc", "TACGtacg"))


def convertRNA(seq) -> str:
    """
    Return the conversion of a RNA sequence to DNA

    :param str: seq: RNA sequence to convert
    :rtype: str
    """

    return seq.translate(str.maketrans("Uu", "Tt"))


def complement(seq) -> str:
    """
    Return the complement of a sequence

    :param str: seq: Sequence to complement
    :rtype: str
    """

    return seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))


def look_forward(iterable, start: int, char: str) -> int:
    """
    Look ahead in an iterable for the next point where a character is different

    :param iterable: iterable: A list/tuple to look forward through
    :param int: start: The initial index to being from
    :param str: char: The character to look for the end of in the sequence
    """

    idx = start
    end_idx = None
    while not end_idx and idx < len(iterable):
        idx += 1
        if iterable[idx] != char:
            end_idx = idx

    return end_idx
