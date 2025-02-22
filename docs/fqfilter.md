# fastq_filtering

Provided a single, or list of, regular expressions, filter out all reads from a given fastq file which match said expression(s).

One use case could be dropping reads with a large amount of poly-A sequences:

`python fqfilter.py -f Sample01.fastq -d "[A]{25}"`

Or perhaps you have a specific adapter at the start of the read that you'd like to isolate

`python fqfilter.py -f Sample01.fastq -k "^ATCGGCTA"`

And you can also do both at once!

`python fqfilter.py -f Sample01.fastq -k "^ATCGGCTA" -d "[A]{25}"`

## Options

|||
|-|-|
| `-f`/`--fastq`        |   Input fastq file to filter                              |
| `-p`/`--pair`         |   Input pair of fastqs. R1 then R2. Files will be syncronized in the ouput for safety of downstream processing!    |
| `-k`/`--keep`         |   Regex expression(s) to keep                             |
| `-d`/`--drop`         |   Regex expression(s) to drop                             |
| ~~`-s`/`--sync`~~         |   ~~Synchronize R1/R2 fastq's for downstream~~ processing     |
| `-v`/`--verbose`      |   Debug output                                            |

## Requirements

Built on Python 3.11, but may be fine back as far as 3.7  
Needs the [BioPython](https://biopython.org/) package
