# SeqTK Subsampling Tool

Generate random subsamples of FASTQ files.

## Overview

This tool wraps SeqTK's `sample` function to create random subsampled FASTQ files. Useful for reducing sequencing depth for testing, analysis, or cost reduction. Supports both single files and batch processing from a TSV index file.

## Usage

```bash
python subsample.py -f <fastq_file> [-n <number>] [-v]
```

## Arguments

| Argument | Short | Long | Type | Default | Required | Description |
| -------- | ----- | ---- | ---- | ------- | -------- | ----------- |
| Input file | `-f` | `--file` | Path | | Yes | FASTQ file or TSV index file |
| Number of reads | `-n` | `--number` | Int | 2000 | No | Target read count for subsample |
| Verbose | `-v` | `--verbose` | Bool | False | No | Enable verbose status messages |

## Input Format

### Single FASTQ File

```bash

# Uncompressed FASTQ
python subsample.py -f reads.fastq -n 5000

# Gzip-compressed FASTQ
python subsample.py -f reads.fastq.gz -n 5000
```

### Batch TSV Index (In Development)

```bash
python subsample.py -f sample_list.tsv -n 2000
```

**TSV Format:**

```tsv
path/to/reads1.fastq.gz
path/to/reads2.fastq.gz
samples/reads3.fastq

# Comments (starting with #) are ignored
```

Paths are relative to TSV file location.

## Supported File Formats

- `*.fastq` - Uncompressed FASTQ
- `*.fastq.gz` - Gzip-compressed FASTQ
- `*.tsv` - Index file with list of FASTQ files (one per line)

## Output

Creates subsampled FASTQ file named:

`subsample_{original_name}.fastq`

**Example:**

```text
Input:  reads.fastq.gz → Output: subsample_reads.fastq.gz
Input:  sample.fastq → Output: subsample_sample.fastq
```

Output preserves input compression (gzip if input .gz).

## SeqTK Integration

Uses SeqTK's `sample` command internally:

```bash
seqtk sample <input_fastq> <fraction>
```

- Requires SeqTK to be installed and in PATH
- Sampling is pseudorandom (consistent across runs with same seed)
- All quality scores and metadata preserved

## Examples

```bash

# Subsample to 5000 reads (default: 2000)
python subsample.py -f large_sample.fastq.gz -n 5000

# Create small test set (1000 reads)
python subsample.py -f whole_genome.fastq -n 1000 -v

# Batch process multiple files (TSV index)
python subsample.py -f sample_list.tsv -n 10000

# Deep read depth to shallow
python subsample.py -f 100M_reads.fastq.gz -n 50000
```

## Processing Steps

1. **Validate input file exists**
2. **Determine input type** (FASTQ or TSV)
3. **If TSV**: Parse file to get list of FASTQ files
4. **For each FASTQ file**:
   - Generate random seed
   - Calculate sampling fraction (desired_reads / total_reads)
   - Run SeqTK sample command
   - Create output file with prefix

## Error Handling

- Validates input file exists before processing
- Skips files listed in TSV that don't exist (warning logged)
- Raises error for invalid file extensions
- Logs parsing warnings for malformed TSV entries

## Logging

Enable detailed output with `-v` flag:

```bash
python subsample.py -f reads.fastq.gz -n 1000 -v
```

Logs include:

- File validation steps
- Input file type detection
- Parse progress for TSV files
- SeqTK command execution

## Requirements

- SeqTK must be installed: `brew install seqtk` or `apt install seqtk`
- BioPython (for FASTQ validation, optional)

## Notes

- **Random seed**: SeqTK uses pseudorandom sampling (reproducible)
- **Memory efficient**: Streams reads rather than loading entire file
- **TSV parsing**: Currently a TODO/incomplete feature
- **Quality preservation**: All quality scores retained in output
- **Large files**: Performance depends on SeqTK efficiency
- **Typical size**: 2000 reads (∼200 Kb) quick validation set

## Typical Use Cases

- **Quality control**: Generate small test set for QC pipeline
- **Benchmarking**: Test analysis on reduced data
- **Validation**: Verify pipeline before processing full sequencing run
- **Cost reduction**: Analyze existing data at lower depth
- **Troubleshooting**: Isolate issues with manageable data size
