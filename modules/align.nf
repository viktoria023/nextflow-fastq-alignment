process ALIGN {
    tag "${sample}"
    publishDir "${params.outputDir}/alignments", mode: 'copy'

    input:
    path reference
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.aligned.bam"), path("${sample}.aligned.bam.bai"), emit: bam
    path "${sample}.flagstat.txt", emit: stats

    script:
    def preset = params.readType == "short" ? "sr" : "map-ont"
    def read_input = params.readType == "short" ? "${reads[0]} ${reads[1]}" : "${reads}"

    """
    # Align reads using minimap2 using the preset for short or long reads
    minimap2 -ax ${preset} -t ${task.cpus} ${reference} ${read_input} | \\
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