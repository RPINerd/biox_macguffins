# Graph-Based Sequence Assembler

De novo sequence assembly using a graph-based overlap approach with Depth-First Search.

## Overview

This tool implements sequence assembly by modeling reads as vertices in a graph where edges represent overlaps. The algorithm finds an assembly path through the graph where each read is visited exactly once, using Depth-First Search (DFS) for path finding.

## Algorithm

### Graph Construction

1. Reads are represented as vertices (graph nodes)
2. Edges represent overlaps between reads:
   - Directed edge from read A to read B if B's sequence can be glued to A's right end
   - Edge weight = protruding length after gluing

### Overlap Detection

Two reads R1 and R2 can be overlapped if:

- Overlap region ≥ half the length of both reads
- Tail of R1 (second half) matches beginning of R2 (first half)
- Overlap can be partial or complete containment

### Assembly Path Finding

- Uses Depth-First Search (DFS) to find Hamiltonian path
- Each read visited exactly once
- Path represents assembly order
- Final sequence constructed by concatenating reads with appropriate overlaps

## Classes

### `Read`

Represents a single sequence read in the graph.

**Attributes:**

- `overlaps` (dict): Map of reads that can follow this one → protruding length
- `visited` (int): Traversal counter (for DFS)
- `visit_limit` (int): Maximum visits allowed (handles duplicate reads)

### `SequenceAssembler`

Main assembly engine.

**Attributes:**

- `reads` (dict): All reads keyed by sequence string
- `path` (list): Assembly order (list of reads)
- `sequence` (str): Final assembled sequence
- `num_reads` (int): Total reads input

**Key Methods:**

- `add_read(read)` - Add sequence to graph
- `read_fasta(handle)` - Load reads from FASTA file
- `calculate_overlap(r1, r2)` - Compute overlap between two reads
- `find_path()` - DFS path finding (inherited from algorithm)
- `assemble()` - Construct final sequence from path

## Usage

```python
from graph_assembler import SequenceAssembler

# Create assembler
assembler = SequenceAssembler()

# Load reads from FASTA
with open("reads.fasta") as f:
    assembler.read_fasta(f)

# Perform assembly
assembler.assemble()

# Access results
print(assembler.sequence)      # Final assembled sequence
print(assembler.path)          # Assembly order
```

## Input Format

Standard FASTA format:

```fasta
>read_1
ACGTACGTACGT...
>read_2
GTACGTACGTAC...
```

## Output

- `sequence`: The final assembled chromosome/contig
- `path`: List of reads in assembly order
- `reads`: Dictionary of all reads with overlap information

## Example

```text
Input reads:
  R1: ATCGGCCAT
  R2: GCCATCGG
  R3: TCGGGCTA

Overlaps detected:
  R1 → R2 (overlap: 5bp, protrusion: +3)
  R2 → R3 (overlap: 4bp, protrusion: +4)

Path found:
  R1 → R2 → R3

Output:
  ATCGGCCATCGGGGCTA
```

## Important Considerations

### Challenges Handled

- **Duplicate reads**: Visit limit tracks identical sequences
- **Branching paths**: DFS explores alternatives until successful path found
- **Contained reads**: Algorithm handles reads completely contained in others

### Limitations

- Greedy overlap matching without quality scoring
- No handling of sequencing errors
- May fail with highly repetitive sequences
- Performance limited by all-vs-all overlap comparison

### When to Use

- Small genomes or contigs (< a few Mb)
- Relatively error-free reads (PacBio HiFi, Illumina long reads)
- Validation/manual assembly review recommended
- Better for targeted assembly than whole-genome

## Notes

- Graph traversal keeps state on each read's visit count
- Protrusion value can be negative (containment scenarios)
- Algorithm is NP-hard; performance degrades with input size
- Consider as experimental/academic implementation
