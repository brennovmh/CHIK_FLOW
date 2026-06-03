#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHIKFLOW } from './workflows/chikflow'

workflow {
    if (params.help) {
        log.info """
        CHIK-FLOW

        Usage:
          nextflow run . --input samplesheet.csv --outdir results -profile docker

        Required arguments:
          --input    CSV samplesheet with columns: sample,fastq_1,fastq_2
          --outdir   Output directory

        Optional arguments:
          --skip_fastqc   Skip FastQC before and after trimming
          --skip_fastp    Skip fastp trimming
          --skip_multiqc  Skip MultiQC aggregation
        """.stripIndent()
        return
    }

    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }

    CHIKFLOW()
}
