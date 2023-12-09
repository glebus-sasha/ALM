#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cpus = 20
params.memory = 40
params.input = '/media/alexandr/KINGSTON/meta_SPAdes/tren.fasta'
params.output = '/home/alexandr/Documents/ALM/data/results'

process GAPSEQ{
    debug true
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/gapseq"

    input:
    path fasta

    output:
    path "*.-all-Pathways.tbl*"
    """
    gapseq find -p all -b 100 -m Bacteria ${fasta}
    """
}
          
workflow{
    GAPSEQ(params.input)
}

