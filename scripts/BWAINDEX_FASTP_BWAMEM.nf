#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cpus = 20
params.memory = 40
params.ref = '/home/oxkolpakova/data/references/genomic.gtf'
params.input = '/home/oxkolpakova/data/raw/*R{1,2}.fq.gz'
params.output = '/home/oxkolpakova/data/results'
params.database = '/srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB'

process FASTP{
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/fastp"

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
    fastp -q 20 -l 140 --trim_poly_g --thread ${task.cpus} \
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
    publishDir "${params.output}/bwaindex"
    
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
    cpus params.cpus
    publishDir "${params.output}/bwamem"
    
    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    
    output:
    path "${sid}.unaligned.bam"
    script:
    """
    bwa mem -t -k 26 ${task.cpus} ${reference} ${reads1} ${reads2} \
    | samtools view -b -f 4 > ${sid}.unaligned.bam 
    """
}

process KRAKEN2 {
    debug true
    cpus params.cpus
    publishDir "${params.output}/kraken2"
    
    input:
    tuple val(sid), path(reads1), path(reads2)
    path database
    
    output:
    path "${sid}kraken2_result.txt"
    path "${sid}kraken2_report.txt"
    
    script:
    """
    kraken2 \
    --db $database \
    --threads ${task.cpus} \
    --output ./${sid}kraken2_result.txt \
    --report ./${sid}kraken2_report.txt \
    --report-zero-counts \
    --use-names \
    --memory-mapping \
    --paired \
    --minimum-base-quality 20 \
    --gzip-compressed \
    ${reads1} ${reads2}
    """
}


if (params.input != false) {
            Channel.fromFilePairs(params.input, checkIfExists: true )
                .set { input_fastqs }
        }
        
workflow{
//    BWAINDEX(params.ref)
    FASTP(input_fastqs)
//    BWAMEM(FASTP.out[0], params.ref, BWAINDEX.out)
    KRAKEN2(FASTP.out[0], params.database)
}

