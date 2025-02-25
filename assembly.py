"""
Problem of assembling the original chromosome (sequence) from
its multiple fragments (reads) is represented with a graph,
where vertices are individual reads and edges are overlaps
between reads. The assembly of the original sequence
is equivalent to finding such a path through the graph
that each read is only used once. The path is found using
the Depth First Search (DFS) algorithm.
"""

from itertools import combinations
from pathlib import Path


class Read:

    """
    Class representing reads in the graph

    Parameters
    ----------
    read : str
        The string that represents a read,
        read from an input file in FASTA format.

    Attributes
    ----------
    overlaps : dict
        The dictionary that holds information about reads that
        can be glued to the current read on the right side.
        Key --- string representing the other read.
        Value --- number of characters the other read protrudes
        with respect to the current read after they are glued.
        For example,

        ATCGGCCAT
            GCCATCGG

        GCCATCGG can be glued to ATCGGCCAT on its right side,
        and protrudes by three characters. The value can be
        negative, in the following example, the value is -2:

        ATCGGCCAT
         TCGGCC

    visited : int
        Number of times this read was visited during
        the graph traversal.

    visit_limit: int
        Limit on the number of times this read can be visited
        during the graph traversal. It is possible by accidence
        to have several reads that are equal to each other,
        but still have a unique way to glue them. The limit is the
        number of reads equal to the current read, including itself.
    """

    def __init__(self, read):
        self.overlaps = {}
        self.visited = 0
        self.visit_limit = 1


class SequenceAssembler:

    """
    Class for sequence assembler

    Attributes
    ----------

    reads : dict
        The dictionary that holds reads.
        Key --- string representing a read.
        Value --- object of Read class.

    path : list
        The list that holds reads in the order they
        should be glued into the origianl sequence.

    sequence : str
        A string representing the original sequence
        to be assembled.

    num_reads : int
        Total number of reads added to the graph.

    """

    def __init__(self):
        self.reads = {}
        self.path = []
        self.sequence = ""
        self.num_reads = 0

    def add_read(self, read):
        """
        Add read to the graph.

        If read is already in the dictinary of reads,
        increment its visit limit.
        """
        r = Read(read)
        if read not in self.reads:
            self.reads[read] = r
        else:
            self.reads[read].visit_limit += 1
        self.num_reads += 1

    def read_fasta(self, handle):
        """
        Read fragments from input file handle.

        For example,

        s = SequenceAssembler()
        with open("data.fasta") as data:
            s.read_fasta(data)

        """
        read = ""
        for line in handle:
            if line[0] == ">":
                if read:
                    self.add_read(read)
                    read = ""
            else:
                read += line.strip()
        self.add_read(read)

    def calculate_overlap(self, r1, r2):
        """
        Check if r1 and r2 can be glued.

        Calculate how much one of them protrudes
        with respect to another after they are glued.
        """
        # We know that reads that can be glued,
        # share at least half of their length.
        # Make sure one is not shorter than
        # the half of another.

        if len(r1) / 2 + len(r1) % 2 <= len(r2) and len(r2) / 2 + len(r2) % 2 <= len(r1):
            # prepare second halves for overlap pre-check
            r1half = int(len(r1) / 2)
            r2half = int(len(r2) / 2)
            tail1 = r1[r1half:]
            tail2 = r2[r2half:]

            # case 1: r1 contains r2 completely
            #
            # For example,
            #
            # ATCGCCGGAT
            #  TCGCCGGA

            pos = r1.find(r2)
            if pos != -1:
                self.reads[r1].overlaps[r2] = pos + len(r2) - len(r1)

            # case 2: r2 contains r1 completely
            #
            # For example,
            #
            #  TCGCCGGA
            # ATCGCCGGAT

            pos = r2.find(r1)
            if pos != -1:
                self.reads[r2].overlaps[r1] = pos + len(r1) - len(r2)

            # case 3: end of r1 overlaps with beginning of r2
            #
            # For example,
            #
            # ATCGCCGGAT
            #  TCGCCGGATGC
            #
            # First check that at least half of r1 is in r2
            # If there is a match, calculate the expected length
            # of overlap and check if they indeed overlap.

            pos = r2.find(tail1)
            if pos != -1:
                overlap = pos + len(tail1)
                if r1[-overlap:] == r2[:overlap]:
                    self.reads[r1].overlaps[r2] = len(r2) - overlap

            # case 4: end of r2 overlaps with beginning of r1
            #
            # For example,
            #
            #   CGCCGGATCC
            #  TCGCCGGAT
            #
            # First check that at least half of r2 is in r1
            # If there is a match, calculate the expected length
            # of overlap and check if they indeed overlap.

            pos = r1.find(tail2)
            if pos != -1:
                overlap = pos + len(tail2)
                if r2[-overlap:] == r1[:overlap]:
                    self.reads[r2].overlaps[r1] = len(r1) - overlap

    def find_path(self, num_visited, read):
        """
        Implements the DFS algorithm.

        For each visited read, we check what reads can be visited next
        and visit one of them. If at some point we are at a dead end,
        we go back and try to visit another one. Continue until we
        visit all reads.
        """
        self.path.append(read)
        r = self.reads[read]
        r.visited += 1
        if num_visited < self.num_reads:
            finished = False
            for other_read in r.overlaps:
                if not finished and self.reads[other_read].visited < self.reads[other_read].visit_limit:
                    finished = self.find_path(num_visited + 1, other_read)
            if not finished:
                self.path.pop()
                r.visited -= 1
        else:
            finished = True
        return finished

    def assemble(self):
        """
        Assemble the original sequence.

        After building the graph (reading all reads)
        calculate all overlaps run the DFS algorithm.
        After finding the path to visit all reads only once,
        assemble the original sequence.
        """
        # Calculate overlaps between each pair of reads.

        for r1, r2 in combinations(self.reads, 2):
            self.calculate_overlap(r1, r2)

        # If there are equal reads, they overlap too

        for read in self.reads:
            if self.reads[read].visit_limit > 1:
                self.reads[read].overlaps[read] = 0

        # Find the read to start the DFS algorithm,
        # The good candidate is a read that can't be glued
        # to any other read on the right side.

        start_candidates = self.reads.copy()

        for read in self.reads:
            r = self.reads[read]
            for other_read in r.overlaps:
                if other_read in start_candidates:
                    del start_candidates[other_read]

        if len(start_candidates):
            for read in start_candidates:
                if len(self.reads[read].overlaps):
                    self.find_path(1, read)
                    break
        else:
            # If there no good candidates where to start
            # the DFS algorithm, try each node.

            for read in self.reads:
                if len(self.reads[read].overlaps):
                    self.find_path(1, read)
                    if len(self.path) == self.num_reads:
                        break

        # Assemble the original sequence:
        # start from the first node in the path,
        # glue subsequent reads, according to how
        # much they are supposed to protrude.

        self.sequence = self.path[0]

        if len(self.path) > 1:
            for i in range(len(self.path) - 1):
                r = self.reads[self.path[i]]
                overlap = r.overlaps[self.path[i + 1]]
                if overlap > 0:
                    self.sequence += self.path[i + 1][-overlap:]
                elif overlap < 0:
                    self.sequence = self.sequence[:overlap]


def main() -> None:
    """Example Implementation"""
    s = SequenceAssembler()

    # Read data
    with Path.open("data.fasta") as data:
        s.read_fasta(data)

    # Run assembly
    s.assemble()

    # Access path and resulting sequence
    print(s.path)
    print(s.sequence)


if __name__ == "__main__":
    main()
