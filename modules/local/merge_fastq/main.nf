process MERGE_FASTQ {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/${meta.id}/fastq", mode: params.publish_dir_mode, pattern: "*.fastq.gz"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_merged_R*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    def paired = meta.single_end ? false : true
    def readNumber = { read ->
        def matcher = read.name =~ /.*(?:^|[._-])R?([12])(?:[._-][0-9]+)?\.f(?:ast)?q\.gz$/
        matcher ? matcher[0][1] : null
    }
    def r1 = reads.findAll { readNumber(it) == '1' }
    def r2 = reads.findAll { readNumber(it) == '2' }
    def single = reads

    if (paired) {
        if (!r1 || !r2) {
            error "Could not identify paired FASTQ files for sample ${meta.id}"
        }
        """
        cat ${r1.join(' ')} > ${meta.id}_merged_R1.fastq.gz
        cat ${r2.join(' ')} > ${meta.id}_merged_R2.fastq.gz

        cat > versions.yml <<-END_VERSIONS
        "${task.process}":
            cat: \$(cat --version | head -n 1 | sed 's/cat (GNU coreutils) //g')
        END_VERSIONS
        """
    } else {
        """
        cat ${single.join(' ')} > ${meta.id}_merged_R1.fastq.gz

        cat > versions.yml <<-END_VERSIONS
        "${task.process}":
            cat: \$(cat --version | head -n 1 | sed 's/cat (GNU coreutils) //g')
        END_VERSIONS
        """
    }
}
