process GENERATE_CONSENSUS {
    tag "${sample}"
    publishDir "${params.outputDir}/consensus", mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path reference

    output:
    tuple val(sample), path("${params.outputPrefix}_${sample}.consensus.fasta"), emit: consensus_fasta

    script:
    """
    # Index reference
    samtools faidx ${reference}
    
    # Compress and index VCF
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz

    # Set VCF file variable
    VCF_FILE="${vcf}.gz"

    # Generate consensus using bcftools
    bcftools consensus -f ${reference} -o ${params.outputPrefix}_${sample}.consensus.fasta \$VCF_FILE
    
    # Update header for clarity
    sed -i '' '1s/.*/>MTB_consensus_${sample}/' ${params.outputPrefix}_${sample}.consensus.fasta
    """

    stub:
    """
    touch ${params.outputPrefix}_${sample}.consensus.fasta
    """
}