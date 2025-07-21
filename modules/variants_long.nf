process CALL_VARIANTS_LONG {
    tag "${sample}"
    publishDir "${params.outputDir}/variants", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.variants.vcf"), emit: vcf

    script:
    """
    # Index the FASTA reference file
    samtools faidx ${reference}

    # Call variants with longshot
    longshot --bam ${bam} --ref ${reference} --out ${sample}.variants.vcf
    """
    
    stub:
    """
    touch ${sample}.variants.vcf
    """
}