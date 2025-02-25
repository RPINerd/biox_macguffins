"""
    Functions related to operations on sequence files (fasta/fastq)

    Operations performed in these functions are limited in scope to the sequences themselves and do not
    take into account the read number of a record (i.e. read1, read2)
"""
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator

    from Bio.SeqRecord import SeqRecord


def dg_count(records: "Generator[SeqRecord, None, None]") -> dict[str, int]:
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


def filter_seq(records: "Generator[SeqRecord, None, None]", regexs: list[re.Pattern] | None, drop: bool = True) -> list["SeqRecord"]:
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

    # logging.info(f"Total Input Records: {total_input}")
    # logging.info(f"Dropped Records: {dropped}")
    # logging.info(f"Saved Records: {len(pruned)}")

    return pruned


def subseq_search(records: "Generator[SeqRecord, None, None]", subseqs: list[str]) -> list[str]:
    """"""
    loc = []
    hit_dict = {}
    for subseq in subseqs:
        hit_list = []
        for record in records:
            regex = re.search(subseq, record)
            if regex:
                hit_list.append(record)
                loc.append(regex.span()[0])

        hit_dict[subseq] = hit_list

        for found in hit_list:
            reads.pop(found)

    if args.reporting:
        logging.info("Generating reports..")
        p_read_file = open("primer_report.txt", "w")
        for key in hit_dict:
            p_read_file.write(str(key) + "|" + str(hit_dict[key]) + "\n")
        p_read_file.close()
        logging.info("Primer Report Written..")
        np_read_file = open("nonprimed.seqs", "w")
        for r in reads:
            np_read_file.write(r + "\n")
        np_read_file.close()
        logging.info("Non-Primed Sequences Written..")
    else:
        logging.info("Skipping report generation..")

    total_nonp = 0
    for seq in reads.keys():
        total_nonp += len(reads[seq])
    avg_pos = sum(loc) / len(loc)
    logging.info("Average subseq start position: " + str(avg_pos))
    logging.info("Unique read sequences without subseq match: " + str(len(reads)))
    logging.info("Total reads without subseq match: " + str(total_nonp))


def uniq(records: "Generator[SeqRecord, None, None]") -> dict[str, list[str]]:
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


def window_shopper(records: "Generator[SeqRecord, None, None]", window_size: int, step: int) -> dict[str, int]:
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
