#!/usr/bin/env nextflow

// path to parameters
params.yml = '/home/ec2-user/ZB4171-Mandolin/env.yml'
params.input = '/home/ec2-user/ZB4171-Mandolin/samplesheet.csv'
params.ref = '/home/ec2-user/ZB4171-Mandolin/reference/KU899140.1_sequence.fasta'
params.bowtie2_script = '/home/ec2-user/ZB4171-Mandolin/src/bowtie2_wrapper.sh'
params.lofreq_script = '/home/ec2-user/ZB4171-Mandolin/src/lofreq_wrapper.sh'

// initiate parameters to channels
Channel.fromPath(file(params.input))
        .splitCsv(sep:',')
        .map {row ->
                def sample_ID = row[0]
                def read1 = file(row[1])
                def read2 = file(row[2])
                def strandedness = row[3]

                return [sample_ID,read1,read2]
        }.into{ ch_sample_fastq; ch_sample_fastq2 }

ch_ref = Channel.value(params.ref)

process bowtie2 {
    publishDir path: "${params.publishdir}/5_bowtie2"
    conda "${params.yml}"

    input:
    tuple val(sample_id), path(read1_paired), path(read2_paired) from ch_sample_fastq

    output:
    tuple val(sample_id), path("*.sam") into ch_bowtie2

    script:
    """
    bash ${params.bowtie2_script} ${sample_id} ${read1_paired} ${read2_paired}
    """
}

process lofreq_preprocess {
    publishDir path: "${params.publishdir}/6_lofreq_preprocess"
    conda "${params.yml}"

    input:
    tuple val(sample_id), path(sam) from ch_bowtie2

    output:
    tuple val(sample_id), path("*_accepted.bam") into ch_lofreq_preprocess

    script:
    """
    samtools view -S -b ${sam} | \
    samtools sort -o ${sample_id}_sorted.bam -
    samtools view -bq 1 ${sample_id}_sorted.bam > ${sample_id}_accepted.bam
    samtools index ${sample_id}_accepted.bam
    """
}

process lofreq_vc {
    publishDir path: "${params.publishdir}/7_lofreq_vc"
    conda "${params.yml}"

    input:
    tuple val(sample_id), path(bam) from ch_lofreq_preprocess

    output:
    tuple val(sample_id), path("*_filtered.vcf") into ch_lofreq

    script:
    """
    bash ${params.lofreq_script} ${sample_id} ${bam}
    """
}