#!/usr/bin/env nextflow

// path to parameters
params.yml = '/home/ec2-user/ZB4171-Mandolin/env.yml'
params.input = '/home/ec2-user/ZB4171-Mandolin/samplesheet.csv'
params.ref = '/home/ec2-user/ZB4171-Mandolin/reference/KU899140.1_sequence.fasta'
params.bowtie2_script = '/home/ec2-user/ZB4171-Mandolin/src/bowtie2_wrapper.sh'

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
    tuple val(sample_id), path("*_accepted.bam") into ch_bowtie2

    script:
    """
    bash ${params.bowtie2_script} ${sample_id} ${read1_paired} ${read2_paired}
    """
}

process lofreq {
    publishDir path: "${params.publishdir}/6_lofreq"
    conda "${params.yml}"

    input:
    tuple val(sample_id), path(bam) from ch_bowtie2

    output:
    tuple val(sample_id), path("*_filtered.vcf") into ch_lofreq

    script:
    """
    lofreq indelqual --dindel --ref ${params.ref} --out ${sample_id}_indelqual.bam ${bam}
    samtools index ${sample_id}_indelqual.bam
    lofreq call -f ${params.ref} -o ${sample_id}.vcf -N -B -q 30 -Q 30 --call-indels --no-default-filter ${sample_id}_indelqual.bam
    lofreq filter --no-defaults --snvqual-thresh 70 --indelqual-thresh 70 --sb-incl-indels -B 60 -a 0.2 -i ${sample_id}.vcf -o ${sample_id}_filtered.vcf
    """
}