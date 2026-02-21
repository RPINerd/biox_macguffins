# BioX MacGuffins

![forthebadge](https://forthebadge.com/badges/made-with-crayons.svg)
![forthebadge](https://forthebadge.com/badges/built-with-grammas-recipe.svg)
![forthebadge](https://forthebadge.com/images/badges/60-percent-of-the-time-works-every-time.svg)

A comprehensive collection of bioinformatics tools and scripts for sequence analysis, quality control, and data processing.

**Requirements:** Python 3.13+, BioPython 1.86+

## Core Modules

### Sequence Operations

- **[fastx.py](fastx.py)** - FASTA/FASTQ file operations (sequence-level operations)
  - Degenerate base counting
  - Sequence validation and filtering

- **[readsets.py](readsets.py)** - Read pair operations (R1/R2 aware)
  - Complexity filtering
  - Paired-end sequence analysis

- **[sambams.py](sambams.py)** - SAM/BAM file manipulation
  - CIGAR string filtering

- **[fastx_utils.py](macguffin_utils.py)** - Utility functions
  - Sequence complement/reverse complement (DNA, RNA, IUPAC)
  - FASTA/FASTQ parsing helpers
  - File collection utilities

### Data Structures

- **[macguffin_classes.py](macguffin_classes.py)** - Core classes
  - `Primer` - Genomic primer representation (chromosome, position, coordinates)
  - `RunCollection` - Pipeline run organization (samples, read sets)
  - `RunSet` - Paired-end read set management

- **[configs.py](configs.py)** - Configuration settings

## Tools

### Sequence Analysis

- **[tools/blastools.py](tools/blastools.py)** - BLAST result manipulation
  - Self-hit removal from BLAST reports

- **[tools/ref_filler.py](tools/ref_filler.py)** - Reference sequence expansion
  - Fills alignment gaps by extracting variants from aligned sequences
  - Generates subsequence FASTA files for missing regions

- **[fasta/pairwise_align.py](fasta/pairwise_align.py)** - Pairwise sequence alignment
  - Aligns sequence pairs using Biopython

- **[tools/ann_to_bed.py](tools/ann_to_bed.py)** - Format conversion
  - Converts UCSC RepMask annotation format to BED format

### Data Processing

- **[clstr_splitter.py](clstr_splitter.py)** - CD-HIT cluster parser
  - Splits CD-HIT cluster output into individual cluster files

- **[dircutadapt.py](dircutadapt.py)** - CutAdapt wrapper
  - Batch adapter trimming for FASTQ files in a directory

- **[subsample.py](subsample.py)** - Subsampling tool
  - Generates subsampled FASTQ files using SeqTK
  - Optional bulk processing from TSV file

- **[fastq/extract_index_files.py](fastq/extract_index_files.py)** - Index extraction
  - Extracts index sequences from paired-end FASTQ (generates I1/I2 from R1/R2)
  - Supports gzipped input

- **[download/fasta_subseq_dl.py](download/fasta_subseq_dl.py)** - Sequence download
  - Downloads FASTA subsequences from NCBI using Accession IDs and coordinates

### Miscellaneous

- **[assembly.py](assembly.py)** - De Bruijn-like sequence assembly
  - Graph-based read assembly using greedy overlap matching
  - Implements Depth-First Search path finding

## Archive (Legacy)

The [archive/](archive/) directory contains deprecated/retired tools:

- `fadiff.py` - FASTA difference (superseded by standard tools)
- `fauniq.py` - FASTA unique (superseded by standard tools)
- `fqfilter.py` - FASTQ filtering
- `cigar_filter.py` - CIGAR filtering (see sambams.py)
- `map_accessiontaxid.py` - Accession to TaxID mapping
- `taxid_annotate.py` - TaxID annotation
- `snp_primer_validate.py` - SNP primer validation
- `read_pair_merger.py` - Paired-end merging
- `tabtodb.py` - TAB to database conversion
- `windowshopper.py` - Window-based sequence analysis
- `rpm_prep` - RPM packaging preparation
