# BED Region Expansion Tool

Expand BED regions that fall below a minimum length threshold.

## Overview

This tool takes a BED file and expands any regions that are shorter than a specified minimum length. Expansion is performed symmetrically around the region center (extending equally upstream and downstream).

## Usage

```bash
python bed_expansion.py -i input.bed -n <minimum_length>
```

## Arguments

| Argument | Short | Long | Type | Required | Description |
| --- | --- | --- | --- | --- | --- |
| Input file | `-i` | `--input` | Path | Yes | Input BED file to expand |
| Minimum length | `-n` | `--num` | Int | Yes | Expand regions to this minimum length |

## Input Format

Standard BED format (at minimum 3 columns):

```tsv
chromosome    start    end    [additional_columns]
```

Columns beyond the first 3 are preserved in output.

## Output Format

Creates a new file named `{input_basename}_expanded.bed` with expanded regions.

- Regions already meeting or exceeding the minimum length are unchanged
- Shorter regions are expanded symmetrically
- If symmetric expansion results in < minimum length, endpoint is incremented by 1

## Examples

```bash
# Expand all regions to at least 1000bp
python bed_expansion.py -i regions.bed -n 1000

# Expand to 500bp minimum
python bed_expansion.py -i coordinates.bed -n 500
```

## Input/Output Example

**Input (regions.bed):**

```bed
chr1    5000    5100      feature1
chr1    10000   10200     feature2
chr2    1000    1010      feature3
```

**Output (regions_expanded.bed with -n 500):**

```bed
chr1    4750    5150      feature1
chr1    10000   10200     feature2
chr2    625     1385      feature3
```

In this example:

- Line 1: 100bp region expanded to 500bp (±200bp from center)
- Line 2: 200bp region already ≥ 500bp, unchanged
- Line 3: 10bp region expanded to 500bp

## Notes

- Expansion is always symmetric around the region midpoint
- If upstream expansion would result in negative coordinates, the region may extend into genome start
- Output file is created in the same directory as input file
