process AA_MUTATIONS {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/${meta.id}/variant_calling", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(variants_vcf), path(variants_index), path(reference_fasta)
    path reference_gff

    output:
    tuple val(meta), path("${meta.id}.aa_mutations.csv"), emit: aa_mutations
    path "versions.yml", emit: versions

    script:
    """
    vcf_to_aa_mutations.py \
        --sample-id "${meta.id}" \
        --vcf "$variants_vcf" \
        --reference-fasta "$reference_fasta" \
        --reference-gff "$reference_gff" \
        --output ${meta.id}.aa_mutations.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
