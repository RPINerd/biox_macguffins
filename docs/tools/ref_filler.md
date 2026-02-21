# Reference Sequence Expansion Tool

Fill alignment gaps in reference sequences using variant sequences.

## Overview

This tool processes a multiple sequence alignment (MSA) and identifies gaps in the reference sequence that are filled in one or more variant sequences. It extracts these "patch" regions and creates separate FASTA files for each gap that can be filled — useful for generating probes, primers, or supplementary sequences to complete reference genomes.

## Algorithm

1. **Identify gaps in reference**: Regions with consecutive `-` characters in the reference (first sequence)
2. **Check variants**: See if any other sequences have actual bases at those positions
3. **Extract patches**: For gaps that can be filled, extract the sequence from variants with flanking context
4. **Output files**: Create separate FASTA for each patch

### Visual Example

```tsv
RefSeq      ACCCACG-----TCCATCAA---GAAGTTCGCGTC------GATGATGCAGTGTGCTAGCTCGACACAGC
VarSeq 1    ACCCACGATACGTCCATCAA---GAAGTTCGCGTC------GATGATGCAGTGT------TCGACACAGC
VarSeq 2    ACCCACG-----TCCATCAATTCGAAGTTCGCGTC------GATGATGCAGTGT------TCGACACAGC
VarSeq 3    ACCCACG-----TCCATCAATTCGAAGTTCGCGTCGGTGACGATGATGCAGTGT------TCGACACAGC
                    Patch1     Patch2             Patch3

Output:
>Patch_1 (filled by VarSeq1)
context_upstream + ATACG + context_downstream

>Patch_2 (filled by VarSeq2,3)
context_upstream + TTC + context_downstream

>Patch_3 (filled by VarSeq3)
context_upstream + GGTGAC + context_downstream
```

## Usage

```bash
python ref_filler.py -f alignment.fasta [-g <gap_size>] [-u <upstream>] [-d <downstream>] [-o <output_dir>]
```

## Arguments

| Argument | Short | Long | Type | Default | Description |
| -------- | ----- | ---- | ---- | ------- | ----------- |
| FASTA file | `-f` | `--fasta` | Path | test/ref_filler.fasta | Input MSA file (first sequence = reference) |
| Gap size | `-g` | `--gap` | Int | 150 | Minimum gap size to fill (bp) |
| Upstream | `-u` | `--upstream` | Int | 100 | Bases upstream of gap to include in patch |
| Downstream | `-d` | `--downstream` | Int | 100 | Bases downstream of gap to include in patch |
| Output dir | `-o` | `--output` | Path | ./ | Directory for output patch FASTA files |

## Input Format

Multiple sequence alignment in FASTA format:

```fasta
>Reference_sequence
ACGTACGTACGT-----TACGTACGT---GGGG
>Variant_1
ACGTACGTACGTACGTACGTACGT---GGGG
>Variant_2
ACGTACGTACGT-----TACGTACGTACGGGGG
```

**Requirements:**

- First sequence is treated as reference
- All sequences must be aligned (same length with gaps as `-`)
- Gaps marked with `-` character

## Output Format

Creates one FASTA file per filled gap named `Patch_<N>.fasta` where N is patch number:

```fasta
>Patch_identifier [source_variant] [position_info]
upstream_sequence + filled_sequence + downstream_sequence
```

## Processing Parameters

### Configurable Constants

Edit these variables at the top of the script:

```python
REGION_SIZE = 150          # Minimum gap size to consider
UPSTREAM_ANCHOR = 100      # Bases upstream context
DOWNSTREAM_ANCHOR = 100    # Bases downstream context
```

Also configurable via command-line arguments (override these defaults).

## Examples

```bash

# Basic usage with defaults (150bp+ gaps)
python ref_filler.py -f aligned.fasta

# Fill smaller gaps (100bp minimum)
python ref_filler.py -f aligned.fasta -g 100

# Include more flanking context (±200bp)
python ref_filler.py -f aligned.fasta -u 200 -d 200

# Specify output directory
python ref_filler.py -f aligned.fasta -o ./patches

# Aggressive gap filling with lots of context
python ref_filler.py -f aligned.fasta -g 50 -u 300 -d 300 -o ./all_patches
```

## Example Output

With input alignment and `-u 50 -d 50`:

```fasta
# File: Patch_1.fasta
>Gap_at_position_1234_from_Variant_2
ACGTACGTACGTACGT[NNNN]TACGTACGTACGT

# File: Patch_2.fasta
>Gap_at_position_5678_from_Variant_1
TACGTACGTACGTACGT[GGTGAC]TACGTACGTACGT
```

## Use Cases

- **Genome assembly**: Fill gaps in fragmented reference genomes
- **Probe design**: Create sequences for gap regions from variants
- **Primer design**: Design primers targeting filled regions
- **Sequence validation**: Verify variant-specific sequences
- **Genomic surveying**: Identify what sequences exist in variants

## Important Notes

- Output format assumes patches represent insertions in reference
- Flanking context allows validation that extracted sequences are correct
- Gaps must be ≥ minimum size threshold to be considered (default 150bp)
- Only gaps that can be filled by at least one variant are output
- Position tracking helps with genomic coordinate mapping
- Anchor regions ensure context for primer/probe design

## Dependencies

Requires utility functions from `macguffin_utils`:

- `look_backward_match()` - Find matching regions upstream
- `look_backward_miss()` - Find gaps upstream
- `look_forward_match()` - Find matching regions downstream
- `look_forward_miss()` - Find gaps downstream
