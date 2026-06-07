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
    path "chikscan.alignment.fasta", emit: alignment
    path "chikscan.distance_matrix.csv", emit: distances
    path "chikscan.phylogeny_metadata.csv", emit: metadata
    path "chikscan.tree.nwk", emit: tree
    path "versions.yml", emit: versions

    script:
    """
    build_chikv_phylogeny.py \
        --consensus ${consensuses} \
        --reference "$reference" \
        --alignment chikscan.alignment.fasta \
        --distances chikscan.distance_matrix.csv \
        --metadata chikscan.phylogeny_metadata.csv \
        --tree chikscan.tree.nwk

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
