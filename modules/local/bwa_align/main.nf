process BWA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bwa=0.7.18"
    container "quay.io/biocontainers/bwa:0.7.18--he4a0461_1"

    input:
    tuple val(meta), path(reads), path(reference_fasta)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path "versions.yml", emit: versions

    script:
    def readArgs = reads.join(' ')

    """
    bwa index "$reference_fasta"

    bwa mem \
        -t ${task.cpus} \
        "$reference_fasta" \
        $readArgs \
        > ${meta.id}.sam

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
