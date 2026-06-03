process VALIDATE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet

    output:
    path "validated_samplesheet.csv", emit: validated_samplesheet
    path "versions.yml", emit: versions

    script:
    """
    validate_samplesheet.py \\
        --input "$samplesheet" \\
        --output validated_samplesheet.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
