process CALL_VARIANTS_SHORT {
    tag "${sample}"
    publishDir "${params.outputDir}/variants", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(sample), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.variants.vcf"), emit: vcf

    script:
    """
    # Index the FASTA reference file
    samtools faidx ${reference}

    # Call variants with freebayes
    freebayes -f ${reference} \\
    --ploidy 1 \\
    --min-base-quality 20 \\
    --min-mapping-quality 30 \\
    --min-coverage 10 \\
    ${bam} > ${sample}.variants.vcf
    """

    stub:
    """
    touch ${sample}.variants.vcf
    """
}