process MULTIQC {
    tag "batch_qc"
    label 'process_low'

    conda "bioconda::multiqc=1.21"
    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"

    publishDir "${params.outdir}/batch_qc/multiqc", mode: params.publish_dir_mode

    input:
    path files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    path "versions.yml", emit: versions

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --outdir .

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //g')
    END_VERSIONS
    """
}
