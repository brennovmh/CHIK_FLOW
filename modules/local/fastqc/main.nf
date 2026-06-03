process FASTQC {
    tag "${meta.id}:${stage}"
    label 'process_low'

    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    publishDir "${params.outdir}/${meta.id}/qc/fastqc/${stage}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"

    input:
    tuple val(meta), path(reads)
    val stage

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"), emit: zip
    path "versions.yml", emit: versions

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads.join(' ')}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//g')
    END_VERSIONS
    """
}
