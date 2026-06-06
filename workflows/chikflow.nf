nextflow.enable.dsl = 2

include { VALIDATE_SAMPLESHEET } from '../modules/local/validate_samplesheet'
include { MERGE_FASTQ          } from '../modules/local/merge_fastq'
include { FASTQC as FASTQC_PRE  } from '../modules/local/fastqc'
include { FASTQC as FASTQC_POST } from '../modules/local/fastqc'
include { FASTP               } from '../modules/local/fastp'
include { MULTIQC             } from '../modules/local/multiqc'
include { VALIDATE_REFERENCE_PANEL } from '../modules/local/validate_reference_panel'
include { BWA_ALIGN           } from '../modules/local/bwa_align'
include { SAMTOOLS_BAM_STATS  } from '../modules/local/samtools_bam_stats'
include { SAMTOOLS_DEPTH      } from '../modules/local/samtools_depth'
include { BCFTOOLS_CONSENSUS } from '../modules/local/bcftools_consensus'
include { GENE_COVERAGE      } from '../modules/local/gene_coverage'
include { SAMPLE_SUMMARY     } from '../modules/local/sample_summary'
include { VARIANT_TABLE      } from '../modules/local/variant_table'
include { AA_MUTATIONS       } from '../modules/local/aa_mutations'

workflow CHIKFLOW {
    main:
    ch_versions = Channel.empty()

    VALIDATE_SAMPLESHEET(file(params.input, checkIfExists: true).toAbsolutePath().toString())
    ch_versions = ch_versions.mix(VALIDATE_SAMPLESHEET.out.versions)

    VALIDATE_SAMPLESHEET.out.validated_samplesheet
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                single_end: row.single_end == 'true'
            ]
            def reads = meta.single_end
                ? [file(row.fastq_1, checkIfExists: true)]
                : [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            [meta, reads]
        }
        .set { ch_reads }

    ch_reads
        .groupTuple(by: 0)
        .map { meta, reads ->
            [meta, reads.flatten()]
        }
        .set { ch_grouped_reads }

    MERGE_FASTQ(ch_grouped_reads)
    ch_analysis_reads = MERGE_FASTQ.out.reads
    ch_versions = ch_versions.mix(MERGE_FASTQ.out.versions)

    ch_multiqc_files = Channel.empty()

    if (!params.skip_fastqc) {
        FASTQC_PRE(ch_analysis_reads, 'pre_trim')
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.map { meta, zip -> zip })
        ch_versions = ch_versions.mix(FASTQC_PRE.out.versions)
    }

    if (!params.skip_fastp) {
        FASTP(ch_analysis_reads)
        ch_trimmed_reads = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, json -> json })
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.map { meta, html -> html })
        ch_versions = ch_versions.mix(FASTP.out.versions)
    } else {
        ch_trimmed_reads = ch_analysis_reads
    }

    if (!params.skip_fastqc) {
        FASTQC_POST(ch_trimmed_reads, 'post_trim')
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.map { meta, zip -> zip })
        ch_versions = ch_versions.mix(FASTQC_POST.out.versions)
    }

    def reference_fasta = file(params.reference_fasta, checkIfExists: true).toAbsolutePath().toString()
    def reference_gff = params.reference_gff
        ? file(params.reference_gff, checkIfExists: true).toAbsolutePath().toString()
        : null

    if (!params.skip_reference_prep) {
        VALIDATE_REFERENCE_PANEL(reference_fasta, reference_gff)
        ch_reference_fasta = VALIDATE_REFERENCE_PANEL.out.fasta
        ch_reference_gff = VALIDATE_REFERENCE_PANEL.out.gff
        ch_versions = ch_versions.mix(VALIDATE_REFERENCE_PANEL.out.versions)
    } else {
        ch_reference_fasta = Channel.value(file(reference_fasta, checkIfExists: true))
        if (reference_gff) {
            ch_reference_gff = Channel.value(file(reference_gff, checkIfExists: true))
        }
    }

    if (!params.skip_alignment) {
        BWA_ALIGN(ch_trimmed_reads, ch_reference_fasta)
        ch_versions = ch_versions.mix(BWA_ALIGN.out.versions)

        SAMTOOLS_BAM_STATS(BWA_ALIGN.out.sam)
        ch_versions = ch_versions.mix(SAMTOOLS_BAM_STATS.out.versions)

        if (!params.skip_coverage || !params.skip_consensus) {
            SAMTOOLS_DEPTH(SAMTOOLS_BAM_STATS.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)
        }

        if (!params.skip_consensus) {
            ch_bam_depth = SAMTOOLS_BAM_STATS.out.bam.join(SAMTOOLS_DEPTH.out.depth)
            BCFTOOLS_CONSENSUS(ch_bam_depth, ch_reference_fasta)
            ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

            VARIANT_TABLE(BCFTOOLS_CONSENSUS.out.variants)
            ch_versions = ch_versions.mix(VARIANT_TABLE.out.versions)

            if (reference_gff) {
                AA_MUTATIONS(BCFTOOLS_CONSENSUS.out.variants, ch_reference_fasta, ch_reference_gff)
                ch_versions = ch_versions.mix(AA_MUTATIONS.out.versions)
            }
        }

        if (!params.skip_coverage && reference_gff) {
            GENE_COVERAGE(SAMTOOLS_DEPTH.out.depth, ch_reference_gff)
            ch_versions = ch_versions.mix(GENE_COVERAGE.out.versions)
        }

        if (!params.skip_coverage && !params.skip_consensus && reference_gff) {
            ch_sample_summary_inputs = SAMTOOLS_BAM_STATS.out.flagstat
                .join(SAMTOOLS_BAM_STATS.out.idxstats)
                .join(SAMTOOLS_DEPTH.out.summary)
                .join(GENE_COVERAGE.out.gene_coverage)
                .join(BCFTOOLS_CONSENSUS.out.consensus)
                .join(BCFTOOLS_CONSENSUS.out.low_coverage)
                .join(BCFTOOLS_CONSENSUS.out.variants)

            SAMPLE_SUMMARY(ch_sample_summary_inputs)
            ch_versions = ch_versions.mix(SAMPLE_SUMMARY.out.versions)
        }
    }

    ch_versions
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    if (!params.skip_multiqc) {
        ch_multiqc_input = ch_multiqc_files.mix(ch_collated_versions).collect()
        MULTIQC(ch_multiqc_input)
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    reads = ch_trimmed_reads
    versions = ch_versions
}
