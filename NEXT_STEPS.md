# ChikScan Next Steps

## Current State

The `main` branch is an executable Nextflow DSL2 development pipeline for CHIKV
short-read analysis.

Implemented blocks:

- samplesheet validation and lane merging
- FastQC, fastp, and MultiQC
- reference FASTA/GFF validation
- BWA-MEM alignment
- sorted/indexed BAM and mapping stats
- per-base depth and genome coverage summary
- GFF-derived gene/CDS coverage summary
- masked consensus FASTA with `bcftools`
- nucleotide variant CSV table
- amino-acid mutation CSV table for CDS-overlapping variants
- per-sample summary CSV
- batch-level sample summary CSV
- nearest-reference genotype/lineage assignment
- batch consensus distance matrix and UPGMA Newick tree
- batch HTML and PDF report

## Validation Baseline

Lightweight local checks:

```bash
python3 -m py_compile \
  bin/aggregate_sample_summaries.py \
  bin/assign_chikv_genotype.py \
  bin/build_chikv_phylogeny.py \
  bin/calculate_gene_coverage.py \
  bin/render_chikscan_report.py \
  bin/select_reference.py \
  bin/summarize_sample.py \
  bin/validate_genotype_references.py \
  bin/validate_reference_panel.py \
  bin/validate_samplesheet.py \
  bin/vcf_to_aa_mutations.py \
  bin/vcf_to_table.py

nextflow run . --help
```

Focused Docker smoke test:

```bash
nextflow run . \
  -profile test,docker \
  --outdir /tmp/chikscan-ci-test \
  --skip_fastqc \
  --skip_fastp \
  --skip_multiqc \
  --min_consensus_depth 1
```

Expected core outputs include:

```text
<outdir>/sample_1/bam/sample_1.sorted.bam
<outdir>/sample_1/coverage/sample_1.depth.tsv
<outdir>/sample_1/coverage/sample_1.coverage_summary.csv
<outdir>/sample_1/coverage/sample_1.gene_coverage.csv
<outdir>/sample_1/assembly/sample_1.consensus.fasta
<outdir>/sample_1/variant_calling/sample_1.variants.vcf.gz
<outdir>/sample_1/variant_calling/sample_1.variants.csv
<outdir>/sample_1/variant_calling/sample_1.aa_mutations.csv
<outdir>/sample_1/genotyping/sample_1.genotype.csv
<outdir>/sample_1/reference_selection/sample_1.reference_selection.csv
<outdir>/sample_1/summary/sample_1.summary.csv
<outdir>/batch_reports/sample_summary.csv
<outdir>/batch_reports/phylogeny/chikscan.phylogeny_metadata.csv
<outdir>/batch_reports/phylogeny/chikscan.tree.nwk
<outdir>/batch_reports/chikscan_report.html
<outdir>/batch_reports/chikscan_report.pdf
```

## Recommended Next Implementation

1. Expand and review the curated CHIKV genotype/reference FASTA.
   - The initial panel includes West African, ECSA, IOL, Asian, and 181/25
     vaccine-strain references.
   - Add additional vaccine-like or outbreak references as they become
     operationally relevant.

2. Improve best-reference selection.
   - The current implementation uses exact read k-mer matching.
   - Next improvements should add minimum score thresholds and clearer
     low-confidence reporting.

3. Improve final reports.
   - Add per-sample coverage plots and batch coverage heatmaps.
   - Add richer mutation tables with gene/CDS context.
   - Replace the current SVG summary with a publication-grade rendered tree.

4. Expand automated validation.
   - Keep the focused CI smoke test fast.
   - Add a scheduled or manually triggered full Docker run without skips.
   - Add regression fixtures for genotype assignment and reference selection.

5. Prepare the first release.
   - Pin and document all runtime containers.
   - Finalize the parameter schema and output contract.
   - Run a pilot batch with real CHIKV data and compare against trusted
     external results.
