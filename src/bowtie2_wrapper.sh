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
indexed_esbl="/home/ec2-user/ZB4171-Mandolin/reference/blm_ref"


# Get the last 2 characters of the sample
last_two_chars="${sample: -2}"

# Check the conditions and run bowtie2 accordingly
if [[ $last_two_chars =~ ^(73|74|76|77)$ ]]; then
    echo "Processing sample with last two characters 73, 74, 76, or 77"
    bowtie2 -x "$indexed_esbl" -1 "$fastq_1" -2 "$fastq_2" -S "${sample}_alignment.sam"
    samtools view -S -b "${sample}_alignment.sam" > "${sample}_alignment.bam"
    samtools sort -o "sorted_${sample}_alignment.bam" "${sample}_alignment.bam"
    samtools view -bq 1 "sorted_${sample}_alignment.bam" > "${sample}_accepted.bam"
    samtools index "${sample}_accepted.bam"
elif [ "$last_two_chars" == "72" ]; then
    # Do something else or leave it blank for now
    echo "Processing sample with last two characters 72"
else
    echo "Processing sample with last two characters not 72, 73, 74, 76, or 77"
    bowtie2 -x "$indexed_blm" -1 "$fastq_1" -2 "$fastq_2" -S "${sample}_alignment.sam"
    samtools view -bS "${sample}_alignment.sam" > "${sample}_alignment.bam"
    samtools view -S -b "${sample}_alignment.sam" > "${sample}_alignment.bam"
    samtools view -bq 1 "sorted_${sample}_alignment.bam" > "${sample}_accepted.bam"
    samtools index "${sample}_accepted.bam"
fi

# Add any additional processing steps or output as needed
# Create index for the reference genome (required only once)


