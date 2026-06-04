process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.20"
    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_1"

    publishDir "${params.outdir}/${meta.id}/coverage", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.depth.tsv"), emit: depth
    tuple val(meta), path("${meta.id}.coverage_summary.csv"), emit: summary
    path "versions.yml", emit: versions

    script:
    """
    samtools depth \
        -aa \
        "$bam" \
        > ${meta.id}.depth.tsv

    awk -v sample_id="${meta.id}" '
        BEGIN {
            FS = "\t"
            total_bases = 0
            depth_sum = 0
            covered_1x = 0
            covered_10x = 0
        }
        {
            depth = \$3 + 0
            total_bases++
            depth_sum += depth
            if (depth >= 1) {
                covered_1x++
            }
            if (depth >= 10) {
                covered_10x++
            }
        }
        END {
            mean_depth = total_bases ? depth_sum / total_bases : 0
            breadth_1x = total_bases ? covered_1x / total_bases : 0
            breadth_10x = total_bases ? covered_10x / total_bases : 0

            print "sample_id,total_bases,covered_bases_1x,covered_bases_10x,mean_depth,breadth_1x,breadth_10x"
            printf "%s,%d,%d,%d,%.4f,%.6f,%.6f\\n", sample_id, total_bases, covered_1x, covered_10x, mean_depth, breadth_1x, breadth_10x
        }
    ' ${meta.id}.depth.tsv > ${meta.id}.coverage_summary.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //g')
    END_VERSIONS
    """
}
