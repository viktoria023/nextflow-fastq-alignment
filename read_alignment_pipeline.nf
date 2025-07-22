#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// Include modules
include { ALIGN } from './modules/align.nf'
include { CALL_VARIANTS_LONG } from './modules/variants_long.nf'
include { CALL_VARIANTS_SHORT } from './modules/variants_short.nf'
include { GENERATE_CONSENSUS } from './modules/consensus.nf'

// Define parameters for input and output paths
params.inputFastq = "./fastq_files"
params.outputDir = "./results"
params.referenceFasta = "./reference/H37Rv.fasta"
params.outputPrefix = "MTB_consensus"
params.readType = null

// Validate required parameters
if (!params.readType) {
    error """
    ERROR: --readType parameter is required!
    
    Please specify either:
      --readType short    (for Illumina/short reads)
      --readType long     (for Nanopore/long reads)
    
    Example:
      nextflow run read_alignment_pipeline.nf --readType short --inputFastq ./data --referenceFasta ./ref.fasta --outputPrefix sample1
    """.stripIndent()
}

if (!(params.readType in ['short', 'long'])) {
    error """
    ERROR: Invalid readType '${params.readType}'
    
    Valid options are:
      --readType short
      --readType long
    """.stripIndent()
}


workflow {
    // Channels
    ch_reference = Channel.fromPath(params.referenceFasta).first() // Allows running multiple fastqs with the same reference

    ch_fastq = params.readType == 'short' ? 
    Channel.fromFilePairs("${params.inputFastq}/*_{1,2}.{fastq.gz}", size: 2, flat: false)
           .ifEmpty { error "No paired-end FASTQ files found. Expected pattern: *_1.fastq.gz/*_2.fastq.gz" } :
    Channel.fromPath("${params.inputFastq}/*.fastq.gz")
           .map { file -> tuple(file.simpleName, file) }

    // Align reads to the reference genome using minimap2 for both long and short reads
    ALIGN(ch_reference, ch_fastq)

    // Call variants using longshot for long reads and freebayes for short reads
    ch_variants = params.readType == 'short' ? 
        CALL_VARIANTS_SHORT(ALIGN.out.bam, ch_reference) : 
        CALL_VARIANTS_LONG(ALIGN.out.bam, ch_reference)

    // Generate consensus sequences from the called variants and reference
    GENERATE_CONSENSUS(ch_variants.vcf, ch_reference)
}
