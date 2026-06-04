process BCFTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.20"
    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    publishDir "${params.outdir}/${meta.id}/assembly", mode: params.publish_dir_mode, pattern: "*.consensus.fasta"
    publishDir "${params.outdir}/${meta.id}/assembly", mode: params.publish_dir_mode, pattern: "*.low_coverage.bed"
    publishDir "${params.outdir}/${meta.id}/variant_calling", mode: params.publish_dir_mode, pattern: "*.vcf.gz*"

    input:
    tuple val(meta), path(bam), path(bai), path(depth)
    path reference_fasta

    output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    tuple val(meta), path("${meta.id}.variants.vcf.gz"), path("${meta.id}.variants.vcf.gz.csi"), emit: variants
    tuple val(meta), path("${meta.id}.low_coverage.bed"), emit: low_coverage
    path "versions.yml", emit: versions

    script:
    def minDepth = params.min_consensus_depth

    """
    awk -v min_depth="${minDepth}" 'BEGIN { OFS = "\\t" } \$3 < min_depth { print \$1, \$2 - 1, \$2 }' \
        "$depth" \
        > ${meta.id}.low_coverage.bed

    bcftools mpileup \
        --fasta-ref "$reference_fasta" \
        --output-type u \
        "$bam" \
        | bcftools call \
            --multiallelic-caller \
            --variants-only \
            --output-type u \
        | bcftools norm \
            --fasta-ref "$reference_fasta" \
            --output-type z \
            --output ${meta.id}.variants.vcf.gz

    bcftools index ${meta.id}.variants.vcf.gz

    bcftools consensus \
        --fasta-ref "$reference_fasta" \
        --mask ${meta.id}.low_coverage.bed \
        ${meta.id}.variants.vcf.gz \
        | awk -v sample_id="${meta.id}" 'BEGIN { renamed = 0 } /^>/ && renamed == 0 { print ">" sample_id; renamed = 1; next } { print }' \
        > ${meta.id}.consensus.fasta

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/bcftools //g')
    END_VERSIONS
    """
}
