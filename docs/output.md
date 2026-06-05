# CHIK-FLOW Output

## Current Output Structure

```text
<OUTDIR>/
├── batch_qc/
│   └── multiqc/
├── pipeline_info/
├── reference_panel/
└── <sample_id>/
    ├── assembly/
    ├── fastq/
    │   └── trimmed/
    ├── bam/
    ├── coverage/
    ├── summary/
    ├── variant_calling/
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

### Alignment

- `<sample_id>/bam/*.sorted.bam`
- `<sample_id>/bam/*.sorted.bam.bai`
- `<sample_id>/bam/*.flagstat.txt`
- `<sample_id>/bam/*.idxstats.txt`

### Coverage

- `<sample_id>/coverage/*.depth.tsv`
- `<sample_id>/coverage/*.coverage_summary.csv`
- `<sample_id>/coverage/*.gene_coverage.csv`

### Summary

- `<sample_id>/summary/*.summary.csv`

### Consensus and Variants

- `<sample_id>/assembly/*.consensus.fasta`
- `<sample_id>/assembly/*.low_coverage.bed`
- `<sample_id>/variant_calling/*.variants.vcf.gz`
- `<sample_id>/variant_calling/*.variants.vcf.gz.csi`
- `<sample_id>/variant_calling/*.variants.csv`

### MultiQC

- `batch_qc/multiqc/multiqc_report.html`
- `batch_qc/multiqc/multiqc_report_data/`

### Pipeline Information

- `pipeline_info/execution_report.html`
- `pipeline_info/execution_timeline.html`
- `pipeline_info/execution_trace.txt`
- `pipeline_info/pipeline_dag.html`
- `pipeline_info/software_versions.yml`

### Reference Panel

- `reference_panel/reference.fasta`
- `reference_panel/reference.gff`
- `reference_panel/reference_panel.csv`

## Planned Output Additions

- `<sample_id>/variant_calling/*.aa_mutations.csv`
- `<sample_id>/genotyping/*.genotype.csv`
- `<sample_id>/phylogeny_tree/*.tre`
- `batch_reports/*.html`
- `batch_reports/*.pdf`
