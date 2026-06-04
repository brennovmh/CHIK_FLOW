# CHIK-FLOW

CHIK-FLOW is a Nextflow DSL2 pipeline scaffold for chikungunya virus (CHIKV)
sequencing analysis. The goal is to provide a RSVrecon-level workflow for raw
FASTQ processing, coverage/depth reporting, consensus generation, genotyping,
phylogeny, and PDF/HTML reporting.

This repository currently contains the first executable foundation:

- samplesheet validation
- optional merge of multiple lanes per sample
- raw read FastQC
- fastp trimming/filtering
- post-trim FastQC
- batch MultiQC
- organized per-sample output directories

## Planned Scope

- CHIKV reference selection against a curated reference panel
- alignment to best matched reference
- BAM sorting, indexing, and mapping statistics
- per-base depth and genome/gene coverage reports
- consensus FASTA with configurable masking thresholds
- nucleotide and amino-acid variant reporting
- CHIKV genotype/lineage assignment
- phylogenetic tree generation
- batch and sample reports in CSV, HTML, and PDF

## Input

Create a CSV samplesheet:

```csv
sample,fastq_1,fastq_2
sample_1,/path/sample_1_R1.fastq.gz,/path/sample_1_R2.fastq.gz
sample_2,/path/sample_2_R1.fastq.gz,/path/sample_2_R2.fastq.gz
```

For single-end data, leave `fastq_2` empty.

If a sample has multiple lanes, add one row per lane with the same `sample`
name. CHIK-FLOW will merge them before analysis.

## Run

```bash
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --reference_fasta references/chikv_panel.fasta \
  --reference_gff references/chikv_panel.gff
```

Run a lightweight configuration check:

```bash
nextflow run . --help
```

## Outputs

Current outputs:

```text
<outdir>/
├── batch_qc/
│   └── multiqc/
├── pipeline_info/
├── reference_panel/
└── <sample>/
    ├── fastq/
    │   └── trimmed/
    ├── log/
    │   └── fastp/
    └── qc/
        ├── fastp/
        └── fastqc/
            ├── post_trim/
            └── pre_trim/
```

Reference panel preparation currently validates FASTA records, optional GFF
seqids, and writes:

```text
reference_panel/reference.fasta
reference_panel/reference.gff
reference_panel/reference_panel.csv
```

Future outputs will add:

```text
<sample>/bam
<sample>/assembly
<sample>/coverage
<sample>/variant_calling
<sample>/genotyping
<sample>/phylogeny_tree
batch_reports
```

## Development Status

This is an initial scaffold. The biological CHIKV-specific layers are planned
but intentionally not stubbed as fake analysis. Each module will be added with
test data and verifiable outputs.
