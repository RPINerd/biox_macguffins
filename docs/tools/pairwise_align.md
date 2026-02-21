# Pairwise Sequence Aligner

Align all sequence pairs from two FASTA files and output alignments.

## Overview

This tool performs pairwise alignment of all sequences from two input FASTA files using BioPython's local alignment algorithm. Creates all-vs-all alignments and writes results to a single output file.

## Usage

```bash
python pairwise_align.py <input_file_1> <input_file_2> <output_file>
```

## Arguments

| Position | Argument | Type | Required | Description |
| -------- | -------- | ---- | -------- | ----------- |
| 1 | Input File 1 | Path | Yes | First input FASTA file |
| 2 | Input File 2 | Path | Yes | Second input FASTA file |
| 3 | Output File | Path | Yes | Output alignment file |

## Input Format

Standard FASTA format for both input files:

```fasta
>sequence_id_1
ACGTACGTACGTACGT
>sequence_id_2
ACGTACGTACGTACGT
```

Multiple sequences per file supported.

## Processing

1. **Sequence Loading**:
   - Reads all sequences from both FASTA files
   - Stores sequences and headers

2. **Special R1 Processing**:
   - Sequences ≤ 30bp from input_file_1 are reverse-complemented
   - Longer sequences kept as-is
   - This behavior is hardcoded and application-specific

3. **Alignment**:
   - Local alignment using BioPython PairwiseAligner
   - All-vs-all comparisons (n×m alignments where n and m are sequence counts)
   - Smith-Waterman algorithm applied

## Alignment Parameters

| Parameter | Value | Description |
| --------- | ----- | ----------- |
| Mode | local | Smith-Waterman local alignment |
| Open gap score | -2.0 | Penalty for gap initiation |
| Extend gap score | -0.5 | Penalty per gap extension |
| End gap score (target) | 0.0 | No penalty for gaps at end of subject |
| End gap score (query) | 0.0 | No penalty for gaps at end of query |

## Output Format

One alignment per block with format:

```fasta
>sequence_id_1::sequence_id_2
aligned_query_sequence
aligned_target_sequence

>next_id_1::next_id_2
aligned_query_sequence
aligned_target_sequence
```

Each alignment includes:

- Header with both sequence IDs separated by `::`
- Query aligned sequence
- Target aligned sequence
- Blank line separator

## Examples

```bash

# Basic usage
python pairwise_align.py query.fasta target.fasta output.fasta

# With different naming
python pairwise_align.py primers.fasta templates.fasta primer_alignments.fasta

# All-vs-all alignment
# 3 sequences in file 1 × 2 sequences in file 2 = 6 alignments
python pairwise_align.py file1.fa file2.fa result.fa
```

## Example Output

```fasta
>read_primer::template_1
AC-GTACGTACGT
ACGGTACGTACGT

>read_primer::template_2
ACGTACGTACGTACGT---
ACGTACGTACGTACGT...

>second_primer::template_1
AC-GTACGTAC-----GT
ACGGTACGTACGTACGGT
```

## Dependencies

- BioPython (Bio.Align, Bio.Seq, Bio.SeqIO)

## Notes

- All combinations of sequences aligned (n × m for n and m sequences)
- Output file size scales quadratically with input sequence count
- Short sequences (≤ 30bp) from file 1 reverse-complemented automatically
- Local alignment finds best local match; useful for finding subsequence matches
- Aligned sequences include gaps marked with `-`
- Suitable for:
  - Primer-template interactions
  - Database similarity searching
  - Multi-target alignment validation
