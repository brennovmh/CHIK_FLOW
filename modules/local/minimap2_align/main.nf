process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::minimap2=2.28"
    container "quay.io/biocontainers/minimap2:2.28--he4a0461_3"

    input:
    tuple val(meta), path(reads)
    path reference_fasta

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path "versions.yml", emit: versions

    script:
    def readArgs = reads.join(' ')

    """
    minimap2 \\
        -ax sr \\
        -t ${task.cpus} \\
        "$reference_fasta" \\
        $readArgs \\
        > ${meta.id}.sam

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
