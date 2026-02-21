# CD-HIT Cluster Splitter

Parse and split CD-HIT cluster output into organized files.

## Overview

This tool processes CD-HIT cluster output (.clstr format) and organizes clusters into separate files based on cluster size. Clusters of size 3 are collected together, clusters of size ≤2 are collected separately, and larger clusters get individual files.

## Usage

```bash
python clstr_splitter.py <cluster_file>
```

## Arguments

- **`cluster_file`** (positional, required): Path to CD-HIT .clstr output file

## Input Format

CD-HIT cluster format with genomic region identifiers:

```text
>Cluster 0
0    LENGTH    start...end...    *
1    LENGTH    start...end...    at 90.5%
>Cluster 1
0    LENGTH    chr1:1000-2000    *
1    LENGTH    chr2:5000-6000    at 95.2%
...
```

Expected region format in cluster lines: `>chr.+:[0-9]+-[0-9]+`

## Output Files

Creates the following output files in the same directory as input:

| File | Contents |
| --- | --- |
| `triples.txt` | All regions from clusters of size exactly 3 (one region per line) |
| `remainder.txt` | All regions from clusters of size 1-2 (one region per line) |
| `ClusterN.txt` | Individual files for each cluster with >3 members |

### Output Format

Each line contains a single genomic region identifier:

```text
chr1:1000-2000
chr2:5000-6000
chr3:10000-11000
```

## Examples

```bash
# Process CD-HIT output
python clstr_splitter.py clusters.clstr

# Output structure created:
# - triples.txt (all size-3 clusters combined)
# - remainder.txt (all size-1 and size-2 clusters combined)
# - Cluster_0.txt (first large cluster)
# - Cluster_1.txt (second large cluster)
# etc.
```

## Cluster Organization Logic

1. **Size = 3**: All regions added to `triples.txt`
2. **Size ≤ 2**: All regions added to `remainder.txt`
3. **Size > 3**: Individual file created: `Cluster_N.txt` where N is the cluster ID

## Notes

- Console output shows total number of clusters found
- Prints warnings for lines not matching expected region format
- Cluster IDs in output filenames use the ID from the cluster header
- Output files created in the same directory as input file
- Existing output files will be overwritten
