# Annotation to BED Converter

Convert UCSC RepMask annotation format to BED format with summary descriptions.

## Overview

This tool converts UCSC RepMask annotation files (e.g., from the RepeatMasker database) to a simplified BED format. It filters for specific repeat types and adds a descriptive 4th column with repeat information.

## Usage

```bash
python ann_to_bed.py <input_file>
```

## Arguments

- **`input_file`** (positional, required): Path to the UCSC RepMask annotation file
  - Download from: <https://hgdownload.soe.ucsc.edu/goldenPath/hgXX/database/rmsk.txt.gz>

## Input Format

The input file should contain 17 tab-separated columns from UCSC RepMask:

```text
bin|swScore|milliDiv|milliDel|milliIns|GenoName|genoStart|genoEnd|genoLeft|strand|repName|repClass|repFamily|repStart|repEnd|repLeft|id
```

## Output Format

Creates `out.bed` file with 4 columns:

```tsv
CHROM    START    END    INFO
```

**INFO column format:** `rep_type|total_len|unit_len|rep_len|rep_base`

- `rep_type`: Repeat type (Simple_repeat, Low_complexity, Satellite)
- `total_len`: Total length of repeat region
- `unit_len`: Length of repeat unit (for Simple_repeat), "." for others
- `rep_len`: Number of repeat copies (for Simple_repeat), "." for others
- `rep_base`: Base/motif sequence

## Supported Repeat Types

- **Simple_repeat**: Tandem repeats with calculated unit length and copy count
- **Low_complexity**: Low-complexity sequence regions
- **Satellite**: Satellite DNA regions

Other repeat types are skipped during conversion.

## Examples

```bash

# Basic usage
python ann_to_bed.py rmsk.txt

# With gzipped input (manually decompress first)
gunzip -c rmsk.txt.gz | python ann_to_bed.py /dev/stdin
```

## Output Example

```tsv
chr1    1000    1050    Simple_repeat|50|10|5|ACGTAC
chr1    2000    2100    Satellite|100|.|.|AAAA
chr1    5000    5150    Low_complexity|150|.|.|GGGGGG
```

## Notes

- Output file is always named `out.bed` in the current working directory
- The tool only processes repeat types: Simple_repeat, Low_complexity, and Satellite
- Other repeat types are skipped
