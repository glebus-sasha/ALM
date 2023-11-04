#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cpus = 5
params.memory = 10
params.ref = '/home/alexandr/Documents/next_try_fastp/reference/hg38_first_1000_lines.fa.gz'
params.input = '/home/alexandr/Documents/next_try_fastp/data/*R{1,2}.fq.gz'
params.output = '.'

process FASTP{
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/results/fastp"

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
            file("${sid}.fastp_stats.json")
            file("${sid}.fastp_stats.html")

    script:
    fq_1_trimmed = sid + '_R1_P.fastq.gz'
    fq_2_trimmed = sid + '_R2_P.fastq.gz'
    """
    fastp -q 20 -l 50 --trim_poly_g --thread ${task.cpus} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]}\
    --out1 $fq_1_trimmed \
    --out2 $fq_2_trimmed \
    --json ${sid}.fastp_stats.json \
    --html ${sid}.fastp_stats.html
    """
}

process BWAINDEX{
    memory params.memory
    cpus params.cpus
    publishDir "${params.output}/results/bwaindex"
    
    input:
    path reference

    output:
    path "*.gz*", emit: bwa_idx

    script:
    """
    bwa index $reference
    """
}

process BWAMEM {
    debug true
    cpus 10
    publishDir "${params.output}/results/bwamem"
    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    output:
    path "${sid}.unaligned.bam"
    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads1} ${reads2} \
    | samtools view -b -f 4 > ${sid}.unaligned.bam 
    """
}

if (params.input != false) {
            Channel.fromFilePairs(params.input, checkIfExists: true )
                .set { input_fastqs }
        }
        
workflow{
    BWAINDEX(params.ref)
    FASTP(input_fastqs)
    BWAMEM(FASTP.out[0], params.ref, BWAINDEX.out)
}

