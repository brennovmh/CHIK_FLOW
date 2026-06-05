process SAMPLE_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/${meta.id}/summary", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(flagstat), path(idxstats), path(coverage_summary), path(gene_coverage), path(consensus), path(low_coverage), path(variants_vcf), path(variants_index)

    output:
    tuple val(meta), path("${meta.id}.summary.csv"), emit: summary
    path "versions.yml", emit: versions

    script:
    """
    summarize_sample.py \
        --sample-id "${meta.id}" \
        --flagstat "$flagstat" \
        --idxstats "$idxstats" \
        --coverage-summary "$coverage_summary" \
        --gene-coverage "$gene_coverage" \
        --consensus "$consensus" \
        --low-coverage-bed "$low_coverage" \
        --variants-vcf "$variants_vcf" \
        --output ${meta.id}.summary.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
