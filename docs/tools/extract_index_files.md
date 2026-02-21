# Extract Index Files

Extract Illumina index sequences from paired-end FASTQ files.

## Overview

This tool extracts dual index sequences (I1 and I2) from paired-end FASTQ read identifiers and writes them to separate FASTQ files. Assumes Illumina format where indices are embedded in the read ID as ``:INDEX1+INDEX2``.

## Usage

```bash
python extract_index_files.py -r R1.fastq.gz [-1 I1_output.fastq] [-2 I2_output.fastq]
```

## Arguments

| Argument | Short | Long | Type | Required | Description |
| --- | --- | --- | --- | --- | --- |
| Read file | `-r` | `--read_file` | Path | Yes | Path to R1 FASTQ file (can be .gz) |
| Index 1 output | `-1` | `--index_1` | Path | No | Output path for I1 FASTQ file |
| Index 2 output | `-2` | `--index_2` | Path | No | Output path for I2 FASTQ file |

## Input Format

Expects Illumina format FASTQ read IDs with indices in the last colon-separated segment:

```fastq
@DEVICEID:FLOWCELL:LANE:TILE:X:Y:INDEX1+INDEX2
ACGTACGTACGT...
+
IIIIIIIIIIII...
```

**Example read ID:**

`@A00123:456:H5LV3DRX2:1:1101:1234:1000:ATTGCC+TAACGG`

Where:

- `ATTGCC` = I1 (Index 1)
- `TAACGG` = I2 (Index 2)

## Output Format

Automatic naming (if not specified):

- Input: `sample_R1.fastq.gz` → Output: `sample_I1.fastq.gz` and `sample_I2.fastq.gz`
- R1 in filename is replaced with I1/I2

Or use custom names with `-1` and `-2` flags.

Output format is standard FASTQ:

```fastq
@READ_ID
INDEX_SEQUENCE
+
QUALITY_SCORES
```

**Note:** Quality scores are set to 'I' (ASCII 73, Illumina quality score 40) for all extracted indices.

## Examples

```bash

# Extract indices with automatic naming
python extract_index_files.py -r sample_R1.fastq.gz

# Creates: sample_I1.fastq.gz, sample_I2.fastq.gz

# With gzipped output (automatic)
python extract_index_files.py -r reads/S001_R1.fastq.gz

# Specify custom output paths
python extract_index_files.py -r R1.fastq -1 indices_1.fastq -2 indices_2.fastq

# Mixed: gzipped input with automatic naming
python extract_index_files.py -r paired_R1.fastq.gz
```

## Handling Gzipped Files

Automatically detects and handles gzipped FASTQ files based on `.gz` extension:

```bash
# Input can be gzipped or uncompressed
python extract_index_files.py -r file.fastq.gz     # gzipped
python extract_index_files.py -r file.fastq        # uncompressed
```

## Error Handling

- **Malformed read IDs**: Issues warning and skips records without proper INDEX1+INDEX2 format
- **Output file exists**: Raises error if I1 or I2 output files already exist
- **Duplicate output files**: Raises error if same file specified for both I1 and I2

## Notes

- Script requires read IDs to contain dual indices separated by `+`
- Single-index or different ID formats may fail silently (warnings logged)
- Quality scores for indices are always set to high quality ('I')
- Logging level set to INFO; warnings for malformed records
- Output files created in same directory as input (if auto-naming)
