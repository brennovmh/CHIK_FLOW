<img width="220" alt="ChikScan" src="assets/report/chikscan_logo.png" />

# ChikScan

ChikScan is a Nextflow DSL2 pipeline for chikungunya virus (CHIKV) sequencing
analysis. It processes raw FASTQ files into QC outputs, alignments, coverage
summaries, consensus sequences, variant tables, genotype/source calls,
phylogeny files, and batch HTML/PDF reports.

The pipeline currently includes:

- samplesheet validation
- optional merge of multiple lanes per sample
- raw read FastQC
- fastp trimming/filtering
- post-trim FastQC
- batch MultiQC
- BWA-MEM alignment to validated CHIKV reference
- sorted/indexed BAM with basic mapping statistics
- per-base depth and basic genome coverage summary
- GFF-derived gene/CDS coverage summary
- consensus FASTA with low-depth masking
- nucleotide variant CSV table
- amino-acid mutation CSV table for CDS-overlapping variants
- per-sample CSV summary across mapping, coverage, consensus, and variants
- batch-level sample summary CSV
- nearest-reference CHIKV genotype/lineage CSV
- batch consensus distance matrix and UPGMA Newick tree
- professional batch HTML and PDF reports with coverage plots, source counts,
  alerts, interpretation notes, and phylogeny rendering
- organized per-sample output directories

## Current Focus

- Expand and curate the CHIKV genotype/source reference panel.
- Improve publication-grade phylogeny rendering.
- Add richer per-sample report pages and longitudinal/batch comparison plots.
- Harden automated biological regression tests around known wild/vaccine
  examples.

## Input

Create a CSV samplesheet:

```csv
sample,fastq_1,fastq_2
sample_1,/path/sample_1_R1.fastq.gz,/path/sample_1_R2.fastq.gz
sample_2,/path/sample_2_R1.fastq.gz,/path/sample_2_R2.fastq.gz
```

For single-end data, leave `fastq_2` empty.

If a sample has multiple lanes, add one row per lane with the same `sample`
name. ChikScan will merge them before analysis.

## Run

```bash
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

ChikScan ships with default reference assets:

- `assets/reference/NC_004162.2.fasta`
- `assets/reference/NC_004162.2.gff`
- `assets/reference/chikv_genotype_references.fasta`

You can override them when needed:

```bash
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --reference_fasta references/chikv_panel.fasta \
  --reference_gff references/chikv_panel.gff \
  --genotype_references references/chikv_genotypes.fasta
```

For surveillance-grade genotype and wild/vaccine calls, use a curated
multi-record genotype FASTA whose headers include labels such as
`|genotype=ECSA|lineage=IOL|source=wild` or
`|genotype=Asian|lineage=vaccine-strain|source=vaccine`.
An initial curated panel is provided at
`assets/reference/chikv_genotype_references.fasta`; provenance and
interpretation notes are documented in `docs/reference_panel.md`.

Run a lightweight configuration check:

```bash
nextflow run . --help
```

Development checks and the focused Docker smoke test are documented in
`docs/development.md`.

## Outputs

Current outputs:

```text
<outdir>/
├── batch_qc/
│   └── multiqc/
├── pipeline_info/
├── reference_panel/
├── batch_reports/
└── <sample>/
    ├── fastq/
    │   └── trimmed/
    ├── bam/
    ├── assembly/
    ├── genotyping/
    ├── reference_selection/
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

Reference panel preparation validates FASTA records, optional GFF
seqids, and writes:

```text
reference_panel/reference.fasta
reference_panel/reference.gff
reference_panel/reference_panel.csv
reference_panel/genotype_references.fasta
reference_panel/genotype_reference_panel.csv
```

Batch reports include:

```text
batch_reports/sample_summary.csv
batch_reports/chikscan_report.html
batch_reports/chikscan_report.pdf
batch_reports/chikscan_phylogeny.svg
batch_reports/phylogeny/chikscan.alignment.fasta
batch_reports/phylogeny/chikscan.distance_matrix.csv
batch_reports/phylogeny/chikscan.phylogeny_metadata.csv
batch_reports/phylogeny/chikscan.tree.nwk
```

## Development Status

ChikScan is an executable development pipeline. Genotype and wild/vaccine
source assignment are implemented as nearest-reference comparisons and are only
as complete as the supplied curated genotype FASTA.
