# CHIK-FLOW Output

## Current Output Structure

```text
<OUTDIR>/
├── batch_qc/
│   └── multiqc/
├── pipeline_info/
├── reference_panel/
├── batch_reports/
│   └── phylogeny/
└── <sample_id>/
    ├── assembly/
    ├── fastq/
    │   └── trimmed/
    ├── bam/
    ├── coverage/
    ├── genotyping/
    ├── reference_selection/
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
- `batch_reports/sample_summary.csv`
- `batch_reports/chikflow_report.html`
- `batch_reports/chikflow_report.pdf`
- `batch_reports/chikflow_phylogeny.svg`

The HTML report includes batch wild/vaccine source counts, alert tables,
genotyping status/source tables, gene coverage plots, and a colored phylogeny.
The PDF contains a compact text summary of the same key sections.

### Genotyping

- `<sample_id>/genotyping/*.genotype.csv`

The genotyping CSV reports the nearest reference, genotype, lineage, identity,
wild/vaccine source, distance, comparable bases, ambiguous bases, status, and
notes. If `--genotype_references` is not provided, the pipeline falls back to
`--reference_fasta`; assignments from a single reference are marked in the note
field as nearest-reference only. Wild/vaccine calls require reference FASTA
headers with `source=wild` or `source=vaccine`; otherwise the source is reported
as `unknown`.

### Reference Selection

- `<sample_id>/reference_selection/*.reference_selection.csv`

The reference selection CSV reports k-mer scores for every candidate reference
record and marks the selected record with `selected=true`. The selected FASTA
record is used for alignment and consensus generation.

### Phylogeny

- `batch_reports/phylogeny/chikflow.alignment.fasta`
- `batch_reports/phylogeny/chikflow.distance_matrix.csv`
- `batch_reports/phylogeny/chikflow.phylogeny_metadata.csv`
- `batch_reports/phylogeny/chikflow.tree.nwk`

The phylogeny metadata CSV reports each tree label, whether the record is a
reference or sample, genotype, lineage, wild/vaccine source, nearest reference,
distance, and comparable bases. Tree labels include the inferred source, for
example `sample_1__source-wild` or `sample_1__source-vaccine`.

### Consensus and Variants

- `<sample_id>/assembly/*.consensus.fasta`
- `<sample_id>/assembly/*.low_coverage.bed`
- `<sample_id>/variant_calling/*.variants.vcf.gz`
- `<sample_id>/variant_calling/*.variants.vcf.gz.csi`
- `<sample_id>/variant_calling/*.variants.csv`
- `<sample_id>/variant_calling/*.aa_mutations.csv`

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
- `reference_panel/genotype_references.fasta`
- `reference_panel/genotype_reference_panel.csv`

## Planned Output Additions

- richer report plots
- publication-grade tree rendering
