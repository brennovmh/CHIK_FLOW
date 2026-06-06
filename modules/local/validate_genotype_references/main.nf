process VALIDATE_GENOTYPE_REFERENCES {
    tag "genotype_references"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/reference_panel", mode: params.publish_dir_mode

    input:
    path references

    output:
    path "genotype_references.fasta", emit: fasta
    path "genotype_reference_panel.csv", emit: panel
    path "versions.yml", emit: versions

    script:
    """
    validate_genotype_references.py \
        --fasta "$references" \
        --output genotype_reference_panel.csv

    cp "$references" genotype_references.fasta

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
