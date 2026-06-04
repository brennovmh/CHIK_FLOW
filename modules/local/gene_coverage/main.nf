process GENE_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/${meta.id}/coverage", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(depth)
    path reference_gff

    output:
    tuple val(meta), path("${meta.id}.gene_coverage.csv"), emit: gene_coverage
    path "versions.yml", emit: versions

    script:
    """
    calculate_gene_coverage.py \
        --sample-id "${meta.id}" \
        --depth "$depth" \
        --gff "$reference_gff" \
        --output ${meta.id}.gene_coverage.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
