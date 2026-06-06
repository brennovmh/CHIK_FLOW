process PHYLOGENY {
    tag "batch_phylogeny"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/batch_reports/phylogeny", mode: params.publish_dir_mode

    input:
    path consensuses
    path reference

    output:
    path "chikflow.alignment.fasta", emit: alignment
    path "chikflow.distance_matrix.csv", emit: distances
    path "chikflow.phylogeny_metadata.csv", emit: metadata
    path "chikflow.tree.nwk", emit: tree
    path "versions.yml", emit: versions

    script:
    """
    build_chikv_phylogeny.py \
        --consensus ${consensuses} \
        --reference "$reference" \
        --alignment chikflow.alignment.fasta \
        --distances chikflow.distance_matrix.csv \
        --metadata chikflow.phylogeny_metadata.csv \
        --tree chikflow.tree.nwk

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
