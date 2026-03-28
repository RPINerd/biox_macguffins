# Fusion Read Detector

K-mer based detection of gene fusion reads in paired-end FASTQ files.

## Overview

This tool performs brute-force detection of fusion reads by building k-mer indices from two reference FASTA sets (fusionA and fusionB), removing shared k-mers to reduce false positives from homologous sequence, then scanning every read for fusion evidence.

Two types of evidence are detected:

- **Spanning reads** — a single read contains k-mers from both fusion partners, suggesting the read crosses the fusion junction
- **Chimeric pairs** — one mate maps to fusionA and the other to fusionB, suggesting the fragment bridges the junction without either read individually spanning it

## Usage

```bash
uv run detect_fusions.py \
  --fusionA first_partner.fasta \
  --fusionB second_partner.fasta \
  --fastq-dir fastqs/
```

## Arguments

| Argument | Short | Long | Type | Default | Required | Description |
| -------- | ----- | ---- | ---- | ------- | -------- | ----------- |
| FusionA reference | `-a` | `--fusionA` | Path(s) | | Yes | One or more FASTA files for the first fusion partner |
| FusionB reference | `-b` | `--fusionB` | Path(s) | | Yes | One or more FASTA files for the second fusion partner |
| FASTQ directory | `-i` | `--fastq-dir` | Path | | Yes | Directory containing `*_R1.fastq.gz` and `*_R2.fastq.gz` files |
| K-mer size | `-k` | `--kmer-size` | Int | 20 | No | K-mer length used for matching |
| Max mismatches | | `--max-mismatches` | Int (0 or 1) | 0 | No | Maximum mismatches allowed per k-mer match |
| Output directory | `-o` | `--out-dir` | Path | `<fastq-dir>/../fusion_hits/` | No | Directory to write hit reads as gzipped FASTQs |

## Algorithm

1. Build k-mer sets (forward + reverse complement) from each fusion partner reference
2. Remove k-mers present in both sets to eliminate homologous cross-hits
3. For every R1/R2 pair in the input directory:
   - Scan each read for fusionA k-mers and fusionB k-mers
   - A single read hitting both → spanning read
   - Mates hitting different partners → chimeric pair
4. Report counts; optionally write hit read pairs to output FASTQs

## Output

Printed summary per sample:

```text
<sample>  total_reads=<N>  spanning=<N>  chimeric=<N>
```

If `--out-dir` is provided, hit reads are written as:

```text
<out-dir>/<sample_stem>_fusion_R1.fastq.gz
<out-dir>/<sample_stem>_fusion_R2.fastq.gz
```

## Requirements

- BioPython

## Notes

- K-mers containing `N` are silently skipped during index construction and scanning
- With `--max-mismatches 1`, a single-base substitution brute-force search is applied (way slower)
- Multiple FASTA files can be supplied to `--fusionA` / `--fusionB` to build a composite reference
