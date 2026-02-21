# CutAdapt Directory Wrapper

Batch adapter trimming for paired-end FASTQ files in a directory.

## Overview

This tool is a wrapper around CutAdapt that automatically processes all paired-end FASTQ files in a directory. It applies adapter trimming to both R1 and R2 reads simultaneously.

**Note:** This script has hardcoded parameters for a specific sequencing project format (MMV2-R format) and may require modification for other use cases.

## Usage

```bash
python dircutadapt.py -d <directory> -a <adapters.fasta>
```

## Arguments

| Argument | Short | Long | Type | Required | Description |
| --- | --- | --- | --- | --- | --- |
| Directory | `-d` | `--directory` | Path | Yes | Directory containing FASTQ files |
| Adapters | `-a` | `--adapters` | Path | Yes | FASTA file with adapter sequences to trim |

## Requirements

- CutAdapt must be installed and accessible in PATH
- Input FASTQ files should follow pattern: `MMV2-R{number}_S{number}_L00{1-4}_R{1,2}(_001.fastq.gz)`

## Input File Format

Expected paired-end file naming:

```text
MMV2-R01_S01_L001_R1_001.fastq.gz
MMV2-R01_S01_L001_R2_001.fastq.gz
MMV2-R02_S02_L001_R1_001.fastq.gz
MMV2-R02_S02_L001_R2_001.fastq.gz
```

## Output Files

For each input pair, creates:

- `trimmed_MMV2-R##_S##_L00#_R1_001.fastq` - Trimmed R1 reads (quality ≥ 25bp)
- `trimmed_MMV2-R##_S##_L00#_R2_001.fastq` - Trimmed R2 reads (quality ≥ 25bp)
- `short_MMV2-R##_S##_L00#_R1_001.fastq` - R1 reads removed (too short)
- `short_MMV2-R##_S##_L00#_R2_001.fastq` - R2 reads removed (too short)
- `reports.tsv` - Trimming statistics per sample

## CutAdapt Parameters

The tool runs CutAdapt with these settings:

```bash
cutadapt -j 12 --report=minimal -n 3 -m 25 \
  -a file:<adapters> -A file:<adapters> \
  --too-short-output short_R1 --too-short-paired-output short_R2 \
  --pair-filter=first \
  -o trimmed_R1 -p trimmed_R2 R1.gz R2.gz
```

| Parameter | Value | Description |
| --- | --- | --- |
| Threads | 12 | Parallel processing |
| Times | 3 | Try up to 3 times to match adapter |
| Min length | 25bp | Reads shorter than 25bp discarded |
| Pair filter | first | Filter pair if R1 is too short |

## Examples

```bash

# Trim adapters from all FASTQ files in directory
python dircutadapt.py -d ./fastq_files -a adapters.fasta

# With adapters in different location
python dircutadapt.py -d /data/reads -a /reference/illumina_adapters.fasta
```

## Adapter FASTA Format

```fasta
>Adapter_1
AGATCGGAAGAGC
>Adapter_2
AGATCGGAAGAGCGTCGT
```

## Notes

- Script extracts sample identifiers from filename and logs to `reports.tsv`
- Only processes R1 files; automatically pairs with corresponding R2 file
- Files not matching the expected naming pattern are skipped
- Output files created in the working directory from which script is called
- 12 parallel threads used for CutAdapt (modify in code if needed)
