process VALIDATE_REFERENCE_PANEL {
    tag "reference_panel"
    label 'process_low'

    publishDir "${params.outdir}/reference_panel", mode: params.publish_dir_mode

    input:
    val reference_fasta
    val reference_gff

    output:
    path "reference_panel.csv", emit: panel
    path "reference.fasta", emit: fasta
    path "reference.gff", optional: true, emit: gff
    path "versions.yml", emit: versions

    script:
    def gffArg = reference_gff ? "--gff \"${reference_gff}\"" : ""
    def gffCopy = reference_gff ? "cp \"${reference_gff}\" reference.gff" : "true"

    """
    validate_reference_panel.py \\
        --fasta "$reference_fasta" \\
        $gffArg \\
        --output reference_panel.csv

    cp "$reference_fasta" reference.fasta
    $gffCopy

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
