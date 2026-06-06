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
    path phylogeny_metadata
    path gene_coverages

    output:
    path "chikflow_report.html", emit: html
    path "chikflow_report.pdf", emit: pdf
    path "chikflow_phylogeny.svg", emit: phylogeny_svg
    path "versions.yml", emit: versions

    script:
    """
    render_chikflow_report.py \
        --sample-summary "$sample_summary" \
        --genotypes ${genotypes} \
        --tree "$tree" \
        --phylogeny-metadata "$phylogeny_metadata" \
        --gene-coverages ${gene_coverages} \
        --html chikflow_report.html \
        --pdf chikflow_report.pdf \
        --phylogeny-svg chikflow_phylogeny.svg

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
