process FINAL_REPORT {
    tag "batch_report"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/batch_reports", mode: params.publish_dir_mode

    input:
    path sample_summary
    path genotypes
    path tree

    output:
    path "chikflow_report.html", emit: html
    path "chikflow_report.pdf", emit: pdf
    path "versions.yml", emit: versions

    script:
    """
    render_chikflow_report.py \
        --sample-summary "$sample_summary" \
        --genotypes ${genotypes} \
        --tree "$tree" \
        --html chikflow_report.html \
        --pdf chikflow_report.pdf

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
