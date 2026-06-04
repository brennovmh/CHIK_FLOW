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

        Reference arguments:
          --reference_fasta  FASTA reference panel used by downstream CHIKV modules
          --reference_gff    Optional GFF annotation for reference panel records

        Optional arguments:
          --skip_fastqc   Skip FastQC before and after trimming
          --skip_fastp    Skip fastp trimming
          --skip_multiqc  Skip MultiQC aggregation
          --skip_reference_prep  Skip reference panel validation
          --skip_alignment  Skip BWA-MEM alignment and BAM processing
          --skip_coverage   Skip per-base depth and coverage summary
          --skip_consensus  Skip variant calling and consensus FASTA generation
          --min_consensus_depth  Minimum depth for unmasked consensus bases
        """.stripIndent()
        return
    }

    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }
    if (!params.skip_reference_prep && !params.reference_fasta) {
        error "Missing required parameter: --reference_fasta"
    }

    CHIKFLOW()
}
