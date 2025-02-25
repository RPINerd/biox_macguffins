""""""
from collections.abc import Generator

from Bio.SeqRecord import SeqRecord


def filter_complexity(
    read1: Generator[SeqRecord, None, None],
    read2: Generator[SeqRecord, None, None]) -> dict[str, tuple[SeqRecord, SeqRecord]]:
    """
    Perform a filtering pass on a set of reads to remove low complexity reads

    Args:
        read1 (Generator[SeqRecord]): The generator object from SeqIO.parse(R1.fastq)
        read2 (Generator[SeqRecord]): The generator object from SeqIO.parse(R2.fastq)

    Returns:
        complex_reads (dict[str, SeqRecord]): Dictionary of complex reads, keyed by read ID

    Raises:
        TypeError: If the read IDs do not match between read1 and read2
    """

    def _complexity_filter(record: SeqRecord) -> bool:
        """
        Filter out reads with low complexity

        Args:
            record (SeqRecord): The read to check

        Returns:
            bool: True if the read is low-complexity
        """
        seq_len = len(record.seq)
        if seq_len <= 50:
            return True
        base_counts = {base: record.seq.count(base) for base in ["A", "T", "C", "G", "N"]}
        for _, count in base_counts.items():
            if count >= (seq_len * 0.5):
                return True
        return (base_counts["C"] + base_counts["G"]) >= (seq_len * 0.75)

    complex_reads: dict[str, tuple[SeqRecord, SeqRecord]] = {}
    dropped_reads: int = 0
    for record_r1, record_r2 in zip(read1, read2):
        generic_read_id = record_r1.id.split()[0]
        if record_r2.id.split()[0] != generic_read_id:
            raise TypeError(f"Read ID mismatch: {generic_read_id} != {record_r2.id.split()[0]}")

        if _complexity_filter(record_r1) or _complexity_filter(record_r2):
            dropped_reads += 1
        else:
            complex_reads[generic_read_id] = (record_r1, record_r2)

    print(f"Total Reads: {len(complex_reads) + dropped_reads}")
    print(f"Dropped Reads: {dropped_reads}")
    print(f"Complex Reads: {len(complex_reads)}")

    return complex_reads
