#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// Include modules
include { ALIGN } from './modules/align.nf'
include { CALL_VARIANTS } from './modules/variants.nf'
include { GENERATE_CONSENSUS } from './modules/consensus.nf'

// Define parameters for input and output paths
params.inputFastq = "./fastq_files"
params.outputDir = "./results"
params.referenceFasta = "./reference/H37Rv.fasta"
params.outputPrefix = "MTB_consensus"

workflow {
    // Channels
    ch_reference = Channel.fromPath(params.referenceFasta).first() // Allows running multiple fastqs with the same reference
    ch_fastq = Channel.fromPath("${params.inputFastq}/*.fastq.gz")
                     .map { file -> tuple(file.simpleName, file) }

    // Debug output
    ch_fastq.view { sample, file -> "Processing: ${sample}" }

    // Align reads to the assembled reference using minimap2
    ALIGN(ch_reference, ch_fastq)

    // Call variants with longshot, producing a VCF
    CALL_VARIANTS(ALIGN.out.bam, ch_reference)

    // Generate consensus sequence from the VCF and reference files
    GENERATE_CONSENSUS(CALL_VARIANTS.out.vcf, ch_reference)
}
