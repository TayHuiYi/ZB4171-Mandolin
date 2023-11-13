#!/usr/bin/env nextflow

// path to parameters
params.yml = '/home/ec2-user/ZB4171-Mandolin/env.yml'
params.input = '/home/ec2-user/ZB4171-Mandolin/samplesheet.csv'
params.ref = '/home/ec2-user/ZB4171-Mandolin/reference/KU899140.1_sequence.fasta'

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

process fastqc {
        publishDir path: "${params.publishdir}/1_fastqc", pattern: "*.html"
	container "staphb/fastqc"
        label "process_medium"

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
        label "process_medium"

        input:
        tuple val (sample_id), path (read1), path (read2) from ch_sample_fastq2

        output:
        tuple val (sample_id), path("*_R1_paired.fastq.gz"), path( "*_R2_paired.fastq.gz") into ch_trim_paired
        tuple val (sample_id), path("*_R1_paired.fastq.gz"), path( "*_R2_paired.fastq.gz") into ch_trim_paired2
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

process interleave_fastq {
	publishDir path: "${params.publishdir}/3_interleaved_fastq"
        conda "${params.yml}"
	label "process_high"

	input:
	tuple val(sample_id), path(read1_paired), path(read2_paired) from ch_trim_paired

	output:
	tuple val(sample_id), path("*_interleaved_modified.fastq.gz") into ch_interleaved_fastq

	script:
	"""
	reformat.sh in1=${read1_paired} in2=${read2_paired} out=${sample_id}_interleaved.fastq.gz fastawrap=0
	zcat ${sample_id}_interleaved.fastq.gz | awk '{if(NR%4==1) {split(\$1, a, "/"); \$1 = a[1] (NR%8<=3 ? "/1" : "/2");} print}' | gzip > ${sample_id}_interleaved_modified.fastq.gz 
	"""
}

process mitobim_reconstruction {
	publishDir path: "${params.publishdir}/4_mitobim"
        container "chrishah/mitobim"
	containerOptions = "--user root"
        label "process_high"

	input:
	tuple val(sample_id), path(interleaved_fastq) from ch_interleaved_fastq
	val (ref) from ch_ref

	output:
	tuple val(sample_id), path("SRR*") 

	script:
	"""
	/home/src/scripts/MITObim.pl -start 1 -end 30 -sample ${sample_id} -ref KU899140.1 -readpool ${interleaved_fastq} --quick /newvol/KU899140.1_sequence.fasta --paired --clean > ${sample_id}_log
	mkdir ${sample_id}
	pattern="iteration*"
	output_dir="${sample_id}/"
        mv "${sample_id}_log" \${output_dir}
	mv \${pattern} \${output_dir}
	"""	
}
