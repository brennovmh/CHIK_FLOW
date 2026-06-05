process VARIANT_TABLE {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/${meta.id}/variant_calling", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(variants_vcf), path(variants_index)

    output:
    tuple val(meta), path("${meta.id}.variants.csv"), emit: variants_csv
    path "versions.yml", emit: versions

    script:
    """
    vcf_to_table.py \
        --sample-id "${meta.id}" \
        --vcf "$variants_vcf" \
        --output ${meta.id}.variants.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
