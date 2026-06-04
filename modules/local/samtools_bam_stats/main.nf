process SAMTOOLS_BAM_STATS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.20"
    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_1"

    publishDir "${params.outdir}/${meta.id}/bam", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat
    tuple val(meta), path("${meta.id}.idxstats.txt"), emit: idxstats
    path "versions.yml", emit: versions

    script:
    """
    samtools sort \\
        --threads ${task.cpus} \\
        -o ${meta.id}.sorted.bam \\
        "$sam"

    samtools index \\
        -@ ${task.cpus} \\
        ${meta.id}.sorted.bam

    samtools flagstat \\
        -@ ${task.cpus} \\
        ${meta.id}.sorted.bam \\
        > ${meta.id}.flagstat.txt

    samtools idxstats \\
        ${meta.id}.sorted.bam \\
        > ${meta.id}.idxstats.txt

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //g')
    END_VERSIONS
    """
}
