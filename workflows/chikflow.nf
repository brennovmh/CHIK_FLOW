nextflow.enable.dsl = 2

include { VALIDATE_SAMPLESHEET } from '../modules/local/validate_samplesheet'
include { MERGE_FASTQ          } from '../modules/local/merge_fastq'
include { FASTQC as FASTQC_PRE  } from '../modules/local/fastqc'
include { FASTQC as FASTQC_POST } from '../modules/local/fastqc'
include { FASTP               } from '../modules/local/fastp'
include { MULTIQC             } from '../modules/local/multiqc'

workflow CHIKFLOW {
    main:
    ch_versions = Channel.empty()

    VALIDATE_SAMPLESHEET(file(params.input, checkIfExists: true))
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
