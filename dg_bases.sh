#!/bin/bash

# This script takes all fasta files (*.fasta or *.fa) in the current directory and inspects each for degenerate bases (N). If any are found, the script will print the name of the file and the number of degenerate bases found

for file in *.fasta *.fa; do
    degenerate_bases=$(grep -o 'N' "$file" | wc -l)
    if [[ $degenerate_bases -gt 0 ]]; then
        echo "$file: $degenerate_bases"
    fi
done