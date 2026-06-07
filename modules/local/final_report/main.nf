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
    path logo

    output:
    path "chikscan_report.html", emit: html
    path "chikscan_report.pdf", emit: pdf
    path "chikscan_phylogeny.svg", emit: phylogeny_svg
    path "versions.yml", emit: versions

    script:
    """
    render_chikscan_report.py \
        --sample-summary "$sample_summary" \
        --genotypes ${genotypes} \
        --tree "$tree" \
        --phylogeny-metadata "$phylogeny_metadata" \
        --gene-coverages ${gene_coverages} \
        --logo "$logo" \
        --html chikscan_report.html \
        --pdf chikscan_report.pdf \
        --phylogeny-svg chikscan_phylogeny.svg

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
