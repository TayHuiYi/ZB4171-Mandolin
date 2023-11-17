#!/bin/bash

# Define an array of values
values=("SRR10277273" "SRR10277274" "SRR10277275" "SRR10277276" "SRR10277277")

# Iterate through the array
for v in "${values[@]}"
do
    # Run the prefetch command
    prefetch $v

    # Run the fastq-dump command with specified options
    fastq-dump $v --split-3 --skip-technical --gzip
done
