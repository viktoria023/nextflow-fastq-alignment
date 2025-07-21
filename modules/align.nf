process ALIGN {
    tag "${sample}"
    cpus 4
    publishDir "${params.outputDir}/alignments", mode: 'copy'

    input:
    path reference
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.aligned.bam"), path("${sample}.aligned.bam.bai"), emit: bam

    script:
    """
    # Align reads using minimap2 and sort
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} | \\
        samtools sort -@ ${task.cpus} -o ${sample}.aligned.bam -

    # Index the BAM file
    samtools index ${sample}.aligned.bam

    # Generate alignment statistics
    samtools flagstat ${sample}.aligned.bam > ${sample}.flagstat.txt
    """

    stub:
    """
    touch ${sample}.aligned.bam
    touch ${sample}.aligned.bam.bai
    touch ${sample}.flagstat.txt
    """
}