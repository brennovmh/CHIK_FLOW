process MERGE_FASTQ {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/${meta.id}/fastq", mode: params.publish_dir_mode, pattern: "*.fastq.gz"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    def paired = meta.single_end ? false : true
    def r1 = reads.findAll { it.name ==~ /.*(_R1|_1|R1).*\.f(ast)?q\.gz/ }
    def r2 = reads.findAll { it.name ==~ /.*(_R2|_2|R2).*\.f(ast)?q\.gz/ }
    def single = reads

    if (paired) {
        """
        cat ${r1.join(' ')} > ${meta.id}_R1.fastq.gz
        cat ${r2.join(' ')} > ${meta.id}_R2.fastq.gz

        cat > versions.yml <<-END_VERSIONS
        "${task.process}":
            cat: \$(cat --version | head -n 1 | sed 's/cat (GNU coreutils) //g')
        END_VERSIONS
        """
    } else {
        """
        cat ${single.join(' ')} > ${meta.id}_R1.fastq.gz

        cat > versions.yml <<-END_VERSIONS
        "${task.process}":
            cat: \$(cat --version | head -n 1 | sed 's/cat (GNU coreutils) //g')
        END_VERSIONS
        """
    }
}
