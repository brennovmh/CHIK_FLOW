# CHIK-FLOW Output

## Current Output Structure

```text
<OUTDIR>/
├── batch_qc/
│   └── multiqc/
├── pipeline_info/
└── <sample_id>/
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

## Current Files

### FastQC

- `<sample_id>/qc/fastqc/pre_trim/*_fastqc.html`
- `<sample_id>/qc/fastqc/pre_trim/*_fastqc.zip`
- `<sample_id>/qc/fastqc/post_trim/*_fastqc.html`
- `<sample_id>/qc/fastqc/post_trim/*_fastqc.zip`

### fastp

- `<sample_id>/fastq/trimmed/*.fastp.fastq.gz`
- `<sample_id>/qc/fastp/*.fastp.html`
- `<sample_id>/qc/fastp/*.fastp.json`
- `<sample_id>/log/fastp/*.fastp.log`

### MultiQC

- `batch_qc/multiqc/multiqc_report.html`
- `batch_qc/multiqc/multiqc_report_data/`

### Pipeline Information

- `pipeline_info/execution_report.html`
- `pipeline_info/execution_timeline.html`
- `pipeline_info/execution_trace.txt`
- `pipeline_info/pipeline_dag.html`
- `pipeline_info/software_versions.yml`

## Planned Output Additions

- `<sample_id>/bam/*.bam`
- `<sample_id>/bam/*.bai`
- `<sample_id>/coverage/*.depth.tsv`
- `<sample_id>/coverage/*.gene_coverage.csv`
- `<sample_id>/assembly/*.consensus.fasta`
- `<sample_id>/variant_calling/*.variants.csv`
- `<sample_id>/variant_calling/*.aa_mutations.csv`
- `<sample_id>/genotyping/*.genotype.csv`
- `<sample_id>/phylogeny_tree/*.tre`
- `batch_reports/*.html`
- `batch_reports/*.pdf`
