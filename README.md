# BioX MacGuffins

![forthebadge](https://forthebadge.com/badges/made-with-crayons.svg)
![forthebadge](https://forthebadge.com/badges/built-with-grammas-recipe.svg)
![forthebadge](https://forthebadge.com/images/badges/60-percent-of-the-time-works-every-time.svg)

A collection of bioinformatics tools and scripts for sequence analysis, quality control, and data processing. Built on the back of BioPython to extend functionality and keep downstream scripting cleaner and more efficient.

## Requirements

Python >= 3.13
BioPython >= 1.86

## Modules

### Sequence Operations

- **[fastx.py](fastx.py)**
  - FASTA/FASTQ file operations (sequence-level operations)
  - Degenerate base counting
  - Sequence validation and filtering

- **[readsets.py](readsets.py)**
  - Read pair operations (R1/R2 aware)
  - Complexity filtering
  - Paired-end sequence analysis

- **[sambams.py](sambams.py)**
  - SAM/BAM file manipulation
  - CIGAR string filtering

- **[macguffin_utils.py](macguffin_utils.py)**
  - Utility functions
  - Sequence complement/reverse complement (DNA, RNA, IUPAC)
  - FASTA/FASTQ parsing helpers
  - File collection utilities

### Data Structures

- **[macguffin_classes.py](macguffin_classes.py)**
  - Core classes
  - `Primer` - Genomic primer representation (chromosome, position, coordinates)
  - `RunCollection` - Pipeline run organization (samples, read sets)
  - `RunSet` - Paired-end read set management

- **[configs.py](configs.py)**
  - Configuration settings

### Download and Processing Tools

- **[macguffin_acquire.py](macguffin_acquire.py)**
  - Data acquisition tools
  - NCBI sequence retrieval
  - Ensemble gene download
