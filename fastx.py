"""
    Functions related to operations on sequence files (fasta/fastq)

    Operations performed in these functions are limited in scope to the sequences themselves and do not
    take into account the read number of a record (i.e. read1, read2)
"""
import logging
import re
from collections.abc import Generator
from pathlib import Path

from Bio.SeqRecord import SeqRecord


def dg_count(records: Generator[SeqRecord]) -> dict[str, int]:
    """
    Count the number of degenerate bases in a set of sequences

    Args:
        records (Generator[SeqRecord]): Generator of SeqRecord objects

    Returns:
        dict[str, int]: Dictionary of sequences and the number of degenerate bases found in each sequence

    Raises:
        ValueError: If the input records are not SeqRecord objects
        Exception: If an unexpected error occurs while processing the records
    """
    dg_counts: dict[str, int] = {}
    for record in records:
        try:
            dg_counts[str(record.id)] = sum(1 for base in str(record.seq) if base.upper() not in {"A", "T", "C", "G"})

        except AttributeError:
            raise ValueError(f"The input records must be SeqRecord objects! Found {type(record)}")

        except Exception as e:
            raise Exception(f"An error occurred while processing the records: {e}")

    return dg_counts


def diff(set1: dict[str, str], set2: dict[str, str]) -> dict[str, str]:
    """
    Compare two fastx sets and return the seqs that are not shared between them

    Args:
        set1 (dict[str, str]): Dictionary of sequences from the first set
        set2 (dict[str, str]): Dictionary of sequences from the second set

    Returns:
        dict[str, str]: Dictionary of sequences that are not shared between the two sets
    """
    return {seq: id for seq, id in set1.items() if seq not in set2}


def filter_seq(records: Generator[SeqRecord], regexs: list[re.Pattern], drop: bool = True) -> list[SeqRecord]:
    """
    Filter a set of sequences based on the provided regular expressions

    Args:
        records (Generator[SeqRecord]): Generator of SeqRecord objects
        regexs (list[re.Pattern] | None): List of regular expressions to match sequences against
        drop (bool): Whether to drop or keep sequences that match the regular expressions

    Returns:
        list[SeqRecord]: List of filtered SeqRecord objects

    Raises:
        ValueError: If no patterns are provided
    """
    total_input = 0
    dropped = 0
    pruned = []
    for record in records:
        total_input += 1
        match = any(re.search(reg, str(record.seq)) for reg in regexs)

        # We are dropping matches, or not keeping because seq did not match
        if (drop and match) or (not drop and not match):
            dropped += 1
            continue

        # Seq did not match drop criteria, or matched keep criteria
        if (drop and not match) or (not drop and match):
            pruned.append(record)

        else:
            # ? Can this even be reached
            raise Exception("An unexpected error occurred while filtering sequences")

    logging.info(f"Total Input Records: {total_input}")
    logging.info(f"Dropped Records: {dropped}")
    logging.info(f"Saved Records: {len(pruned)}")

    return pruned


def subseq_search(records: Generator[SeqRecord], subseqs: list[str]) -> None:
    """
    Given a list of query sequences, iterate through a set of records and summarize the presence of the query sequences

    Generates 2 report files:
    - hit_reports.txt: Contains the query sequence and the ids of the records that contain the query sequence
    - nohit_reports.txt: Contains the ids of the records that do not contain any of the query sequences

    Args:
        records (Generator[SeqRecord]): Generator of SeqRecord objects
        subseqs (list[str]): List of subseqs to search for

    Returns:
        None

    Raises:
        ValueError: If no subseqs are provided
        ValueError: If an input record is not a SeqRecord object
    """
    if not subseqs:
        raise ValueError("No subseqs provided!")

    loc = []
    hit_dict: dict[str, list[SeqRecord]] = {}
    nohit_records: list[SeqRecord] = []
    hit_file: Path = Path("hit_reports.txt")
    nohit_file: Path = Path("nohit_reports.txt")
    for record in records:
        if not isinstance(record, SeqRecord):
            raise ValueError(f"Input records must be SeqRecord objects! Found {type(record)}")

        if record.id is None:
            continue  # Skip records without an id

        nohit = True
        for subseq in subseqs:
            regex: re.Match | None = re.search(subseq, record.seq)
            if regex:
                nohit = False
                try:
                    hit_dict[subseq].append(record)
                except KeyError:
                    hit_dict[subseq] = [record]
                loc.append(regex.span()[0])
        if nohit:
            nohit_records.append(record)

    logging.info("Generating reports..")

    with Path.open(hit_file, "w") as hf:
        for query_seq, record_list in hit_dict.items():
            hf.write(f"{query_seq}|{",".join([record.id for record in record_list])}\n")
    logging.info("Primer Report Written..")

    with Path.open(nohit_file, "w") as nhf:
        nhf.write("\n".join([record.id for record in nohit_records]))
    logging.info("Non-Primed Sequences Written..")

    avg_pos = sum(loc) / len(loc)
    logging.info(f"Average subseq start position: {avg_pos}")
    logging.info(f"Reads without subseq hit: {len(nohit_records)}")


def uniq(records: Generator[SeqRecord]) -> dict[str, list[str]]:
    """
    Consolidate a fastx file into dictionary of unique sequences

    Args:
            records (Generator[SeqRecord]): Generator of SeqRecord objects

    Returns:
            dict[str, list[str]]: Dictionary of unique sequences where the key is the sequence and the value is a list of ids

    Raises:
            ValueError: If the input records are not SeqRecord objects
            Exception: If an unexpected error occurs while processing the records
    """
    uniq_seqs: dict[str, list[str]] = {}
    for record in records:
        try:
            uniq_seqs[str(record.seq)].append(str(record.id))
        except KeyError:
            uniq_seqs[str(record.seq)] = [str(record.id)]
        except AttributeError:
            raise ValueError(f"The input records must be SeqRecord objects! Found {type(record)}")
        except Exception as e:
            raise Exception(f"An error occurred while processing the records: {e}")

    return uniq_seqs


def window_shopper(records: Generator[SeqRecord], window_size: int, step: int) -> dict[str, int]:
    """
    Generate a dictionary of sequences and their counts from a sliding window analysis

    Args:
        records (Generator[SeqRecord]): Generator of SeqRecord objects
        window_size (int): Size of the sliding window
        step (int): How many steps to take per window iteration

    Returns:
        dict[str, int]: Dictionary of sequences and their counts
    """
    window_dict: dict[str, int] = {}
    for record in records:
        seq = str(record.seq)
        start = 0
        while start < len(seq) - window_size:
            end = start + window_size
            window = seq[start:end]
            try:
                window_dict[window] += 1
            except KeyError:
                window_dict[window] = 1
            start += step

    return window_dict
