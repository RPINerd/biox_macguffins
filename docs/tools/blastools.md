# BLAST Tools

Utilities for manipulating and processing BLAST output files.

## Overview

This tool provides functions for cleaning and filtering BLAST tabular format results. Currently implements self-hit removal with support for additional filtering modes.

## Usage

```bash
python blastools.py -i input.txt -o output.txt [options]
```

## Arguments

| Argument | Short | Long | Type | Default | Description |
| --- | --- | --- | --- | --- | --- |
| Input file | `-i` | `--input` | Path | Required | Input BLAST tabular format file |
| Output file | `-o` | `--output` | Path | output.txt | Output file path |
| Verbose | `-v` | `--verbose` | Flag | False | Enable verbose output |
| Trim self-hits | | `--trimselfhits` | Flag | False | Remove self-hit matches |
| Pick best hit | | `--pickbesthit` | Flag | False | Keep only best hit per query |

## Input Format

Expects BLAST tabular format with at least 2 columns containing hit coordinates in the format:

```tsv
query_id    subject_id    [other columns...]
```

Coordinate format example: `chr1:1000-2000`

## Functional Modes

### Trim Self-Hits

Removes matches where query and subject are identical (same coordinates).

Uses regex pattern: `^(chr.+:[0-9]+\-[0-9]+)\t(chr.+:[0-9]+\-[0-9]+)`

```bash
python blastools.py -i hits.txt --trimselfhits -o filtered.txt
```

### Pick Best Hit (Incomplete)

Function stub for selecting highest quality match per query. Currently unimplemented.

## Examples

```bash

# Remove self-hits from BLAST results
python blastools.py -i blast_results.txt --trimselfhits -o cleaned.txt

# With verbose output
python blastools.py -i blast_results.txt --trimselfhits -v -o cleaned.txt

# Combine operations
python blastools.py -i blast_results.txt --trimselfhits --pickbesthit -o best_unique.txt
```

## Output

Progress feedback printed to stdout (e.g., "Removing self hits.. X%").

Output file contains all retained hits in original tabular format.

## Notes

- Input file must be in BLAST tabular format (standard output from `blastn -outfmt 6`)
- Coordinate pattern matching is specific to genomic coordinates (chr:start-end format)
- `--pickbesthit` function is not yet implemented
- Progress updates shown as percentage during self-hit removal
