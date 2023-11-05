#!/usr/bin/env nextflow

params.yml = '/home/ec2-user/ZB4171-Mandolin/env.yml'
params.input = '/home/ec2-user/ZB4171-Mandolin/samplesheet.csv'

Channel.fromPath(file(params.input))
        .splitCsv(sep:',')
        .map {row ->
                def sample_ID = row[0]
                def read1 = file(row[1])
                def read2 = file(row[2])
                def strandedness = row[3]

                return [sample_ID,read1,read2]
        }.into{ ch_sample_fastq; ch_sample_fastq2 }

process fastqc {
        publishDir path: "${params.publishdir}/1_fastqc", pattern: "*.html"
	container "staphb/fastqc"
        label "process_high"

        input:
        tuple val (sample_id), path (read1), path (read2) from ch_sample_fastq

        output:
        tuple val (sample_id), path("*.html") 

        script:
        """
        fastqc -t 16 ${read1} ${read2}
        """
}

process trimmomatic {
        publishDir path: "${params.publishdir}/2_trimmomatic"
	conda "${params.yml}"
        label "process_high"

        input:
        tuple val (sample_id), path (read1), path (read2) from ch_sample_fastq2

        output:
        tuple val (sample_id), path("*_R1_paired.fastq.gz"), path( "*_R2_paired.fastq.gz") into ch_trim_paired
        tuple val (sample_id), path("*_R1_unpaired.fastq.gz"), path( "*_R2_unpaired.fastq.gz") into ch_trim_unpaired

        script:
        """
        trimmomatic PE \
        -phred33 \
	-threads 16 \
	-trimlog "${sample_id}_trimlog.txt" \
        ${read1} ${read2} \
        "trimmed_${sample_id}_R1_paired.fastq.gz" "trimmed_${sample_id}_R1_unpaired.fastq.gz" \
        "trimmed_${sample_id}_R2_paired.fastq.gz" "trimmed_${sample_id}_R2_unpaired.fastq.gz" \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 \
        ILLUMINACLIP:/home/ec2-user/ZB4171-Mandolin/TruSeq3-PE.fa:2:30:10 \
        """
}

