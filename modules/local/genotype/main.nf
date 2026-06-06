process GENOTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/${meta.id}/genotyping", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(consensus)
    path references

    output:
    tuple val(meta), path("${meta.id}.genotype.csv"), emit: genotype
    path "versions.yml", emit: versions

    script:
    """
    assign_chikv_genotype.py \
        --sample-id ${meta.id} \
        --consensus "$consensus" \
        --references "$references" \
        --output ${meta.id}.genotype.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
