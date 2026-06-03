process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "quay.io/biocontainers/fastp:0.23.4--hadf994f_3"

    publishDir "${params.outdir}/${meta.id}/fastq/trimmed", mode: params.publish_dir_mode, pattern: "*.fastp.fastq.gz"
    publishDir "${params.outdir}/${meta.id}/qc/fastp", mode: params.publish_dir_mode, pattern: "*.fastp.{html,json}"
    publishDir "${params.outdir}/${meta.id}/log/fastp", mode: params.publish_dir_mode, pattern: "*.fastp.log"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastp.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.html"), emit: html
    tuple val(meta), path("*.fastp.json"), emit: json
    tuple val(meta), path("*.fastp.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = meta.single_end
        ? "-i ${reads[0]} -o ${meta.id}.fastp.fastq.gz"
        : "-i ${reads[0]} -I ${reads[1]} -o ${meta.id}_R1.fastp.fastq.gz -O ${meta.id}_R2.fastp.fastq.gz"

    """
    fastp \\
        --thread ${task.cpus} \\
        $args \\
        --html ${meta.id}.fastp.html \\
        --json ${meta.id}.fastp.json \\
        2> ${meta.id}.fastp.log

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //g')
    END_VERSIONS
    """
}
