#!/usr/bin/env nextflow

params.input = '/home/ec2-user/ZB4171-Mandolin/samplesheet.csv'

Channel.fromPath(file(params.input))
        .splitCsv(sep:',')
        .map {row ->
                def sample_ID = row[0]
                def read1 = file(row[1])
                def read2 = file(row[2])
                def strandedness = row[3]

                return [sample_ID,read1,read2]
        }.set { ch_sample_fastq }

process fastqc {
        publishDir path: "${params.publishdir}/1_fastqc", pattern: "*.html"
	container "staphb/fastqc"
        label "process_low"

        input:
        tuple val (sample_id), path (read1), path (read2) from ch_sample_fastq

        output:
        tuple val (sample_id), path("*.html") 

        script:
        """
        fastqc -t 10 ${read1} ${read2}
        """
}
