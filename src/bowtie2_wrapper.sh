#!/bin/bash

# Check if the correct number of command line arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 sample fastq_1 fastq_2"
    exit 1
fi

# Assign command line arguments to variables
sample="$1"
fastq_1="$2"
fastq_2="$3"

echo "fastq_1: $fastq_1"
echo "fastq_2: $fastq_2"


# Define the reference sequences with absolute paths
indexed_blm="/home/ec2-user/ZB4171-Mandolin/reference/blm_ref"
indexed_esbl="/home/ec2-user/ZB4171-Mandolin/reference/esbl_ref"
indexed_pres="/home/ec2-user/ZB4171-Mandolin/reference/pres_ref"

# Get the last 2 characters of the sample
last_two_chars="${sample: -2}"

# Check the conditions and run bowtie2 accordingly
if [[ $last_two_chars =~ ^(73|74|76|77)$ ]]; then
    echo "Processing sample with last two characters 73, 74, 76, or 77"
    bowtie2 -x "$indexed_esbl" --very-fast --threads 5 -1 "$fastq_1" -2 "$fastq_2" -S "${sample}_alignment.sam"
elif [ "$last_two_chars" == "72" ]; then
    # Do something else or leave it blank for now
    echo "Processing sample with last two characters 72"
    bowtie2 -x "$indexed_pres" --very-fast --threads 5 -1 "$fastq_1" -2 "$fastq_2" -S "${sample}_alignment.sam"
else
    echo "Processing sample with last two characters not 72, 73, 74, 76, or 77"
    bowtie2 -x "$indexed_blm" --very-fast --threads 5 -1 "$fastq_1" -2 "$fastq_2" -S "${sample}_alignment.sam"
fi

# Add any additional processing steps or output as needed
# Create index for the reference genome (required only once)


