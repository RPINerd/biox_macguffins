# Tandem Repeat Inspector

Detect assembled contigs containing tandem repeats of a specified gene.

## Overview

Given a set of assembled contigs (FASTA) and a gene accession, this tool fetches the gene's exon sequences from NCBI and uses approximate string matching to identify contigs that contain multiple exon hits — a signature of tandem gene duplication or repeat expansion.

## Usage

```bash
python detect_tandem_repeats.py -r contigs.fasta -a GENE_ACCESSION [-o output.txt] [-v]
```

## Arguments

| Argument | Short | Long | Type | Required | Description |
| -------- | ----- | ---- | ---- | -------- | ----------- |
| Read file | `-r` | `--read_file` | Path | Yes | Assembled contigs in FASTA format |
| Accession | `-a` | `--accession` | String | Yes | NCBI accession ID of the gene to search for |
| Output | `-o` | `--output` | Path | No | Output file path (default: stdout) |
| Verbose | `-v` | `--verbose` | Flag | No | Enable verbose/debug logging |

## Input Format

Standard FASTA format for assembled contigs:

```fasta
>contig_001
ACGTACGTACGTACGT...
>contig_002
ACGTACGTACGTACGT...
```

## Algorithm

1. Fetch exon sequences for the given accession from NCBI (via `macguffin_acquire.fetch_exons`)
2. For each contig, run approximate alignment of every exon against the contig using `edlib` in HW (infix) mode
3. An exon is considered a hit when alignment similarity ≥ 90% (`MIN_EXON_SIMILARITY = 0.90`)
4. Contigs with one or more exon hits are recorded along with the ordered list of matching exon numbers

## Output

Tab-separated: contig ID followed by a comma-separated list of matched exon numbers.

**To file (`-o output.txt`):**

```tsv
contig_001    exon_1,exon_3,exon_5
contig_007    exon_2,exon_2,exon_4
```

**To stdout (no `-o`):** Same format printed to standard output.

## Requirements

- `edlib` — fast approximate string matching
- BioPython

## Notes

- Exon similarity threshold is hardcoded at 90% (`MIN_EXON_SIMILARITY`)
- A contig appearing multiple times for the same exon number is a strong indicator of a tandem repeat
- Use `-v` to enable per-exon debug output showing match positions and similarity scores
