#!/bin/bash

# Check if the correct number of command line arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 sample bam"
    exit 1
fi

# Assign command line arguments to variables
sample="$1"
bam="$2"

echo "bam: $bam"


# Define the reference sequences with absolute paths
ref_blm="/home/ec2-user/ZB4171-Mandolin/reference/KU899140.1_sequence.fasta"
ref_esbl="/home/ec2-user/ZB4171-Mandolin/reference/esbl_ref.fasta"
ref_pres="/home/ec2-user/ZB4171-Mandolin/reference/pres_ref.fasta"

# Get the last 2 characters of the sample
last_two_chars="${sample: -2}"

# Check the conditions and run bowtie2 accordingly
if [[ $last_two_chars =~ ^(73|74|76|77)$ ]]; then
    echo "Processing sample with last two characters 73, 74, 76, or 77"
    lofreq indelqual --dindel --ref ${ref_esbl} --out ${sample_id}_indelqual.bam ${bam}
    samtools index ${sample_id}_indelqual.bam
    lofreq call -f ${ref_esbl} -o ${sample_id}.vcf -N -B -q 30 -Q 30 --call-indels --no-default-filter ${sample_id}_indelqual.bam
    lofreq filter --no-defaults --snvqual-thresh 70 --indelqual-thresh 70 --sb-incl-indels -B 60 -a 0.2 -i ${sample_id}.vcf -o ${sample_id}_filtered.vcf
elif [ "$last_two_chars" == "72" ]; then
    # Do something else or leave it blank for now
    echo "Processing sample with last two characters 72"
    lofreq indelqual --dindel --ref ${ref_pres} --out ${sample_id}_indelqual.bam ${bam}
    samtools index ${sample_id}_indelqual.bam
    lofreq call -f ${ref_pres} -o ${sample_id}.vcf -N -B -q 30 -Q 30 --call-indels --no-default-filter ${sample_id}_indelqual.bam
    lofreq filter --no-defaults --snvqual-thresh 70 --indelqual-thresh 70 --sb-incl-indels -B 60 -a 0.2 -i ${sample_id}.vcf -o ${sample_id}_filtered.vcf
else
    echo "Processing sample with last two characters not 72, 73, 74, 76, or 77"
    lofreq indelqual --dindel --ref ${ref_blm} --out ${sample_id}_indelqual.bam ${bam}
    samtools index ${sample_id}_indelqual.bam
    lofreq call -f ${ref_blm} -o ${sample_id}.vcf -N -B -q 30 -Q 30 --call-indels --no-default-filter ${sample_id}_indelqual.bam
    lofreq filter --no-defaults --snvqual-thresh 70 --indelqual-thresh 70 --sb-incl-indels -B 60 -a 0.2 -i ${sample_id}.vcf -o ${sample_id}_filtered.vcf
fi
