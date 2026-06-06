# CHIK-FLOW Next Steps

## Current State

The local checkout is synced with GitHub `main`. Use `git log --oneline -5`
for the current commit history.

Latest implemented blocks:

- Real CHIKV RefSeq reference: `NC_004162.2`
- BWA-MEM alignment
- sorted/indexed BAM and mapping stats
- per-base depth and genome coverage summary
- masked consensus FASTA with `bcftools`
- GFF-derived gene/CDS coverage summary
- per-sample summary CSV
- nucleotide variant CSV table
- amino-acid mutation CSV table for CDS-overlapping variants
- batch-level sample summary CSV

The per-sample summary combines mapping stats, genome coverage, GFF feature
coverage highlights, consensus metrics, low-coverage masking, and VCF record
counts.

Outputs:

```text
<sample>/summary/*.summary.csv
<sample>/variant_calling/*.variants.csv
<sample>/variant_calling/*.aa_mutations.csv
batch_reports/sample_summary.csv
```

## Validation Already Completed

Python syntax validation:

```bash
python3 -m py_compile \
  bin/aggregate_sample_summaries.py \
  bin/vcf_to_aa_mutations.py \
  bin/vcf_to_table.py \
  bin/summarize_sample.py \
  bin/calculate_gene_coverage.py \
  bin/validate_reference_panel.py \
  bin/validate_samplesheet.py
```

Nextflow help/config validation:

```bash
nextflow run . --help
```

Focused Singularity test:

```bash
nextflow run . -profile test,singularity \
  --outdir /tmp/chikflow-batch-test \
  --skip_fastqc \
  --skip_fastp \
  --skip_multiqc \
  --min_consensus_depth 1
```

This Singularity run completed successfully and produced:

```text
/tmp/chikflow-batch-test/sample_1/summary/sample_1.summary.csv
/tmp/chikflow-batch-test/sample_1/variant_calling/sample_1.variants.csv
/tmp/chikflow-batch-test/sample_1/variant_calling/sample_1.aa_mutations.csv
/tmp/chikflow-batch-test/batch_reports/sample_summary.csv
```

## Docker Status

Docker is accessible in the current Codex session.

Docker smoke test:

```bash
docker run --rm hello-world
```

Observed result:

```text
Hello from Docker!
```

Focused Docker test:

```bash
nextflow run . -profile test,docker \
  --outdir /tmp/chikflow-batch-docker-test \
  --skip_fastqc \
  --skip_fastp \
  --skip_multiqc \
  --min_consensus_depth 1
```

This Docker run completed successfully and produced:

```text
/tmp/chikflow-batch-docker-test/sample_1/summary/sample_1.summary.csv
/tmp/chikflow-batch-docker-test/sample_1/variant_calling/sample_1.variants.csv
/tmp/chikflow-batch-docker-test/sample_1/variant_calling/sample_1.aa_mutations.csv
/tmp/chikflow-batch-docker-test/batch_reports/sample_summary.csv
```

Note: the first Docker pipeline attempt reached `VARIANT_TABLE` but failed while
pulling `python:3.12` due to a Docker Hub TLS handshake timeout. Re-running the
test succeeded.

For future full Docker checks:

```bash
nextflow run . -profile test,docker --outdir /tmp/chikflow-docker-test
```

Expected Docker test outputs include:

```text
/tmp/chikflow-docker-test/sample_1/bam/sample_1.sorted.bam
/tmp/chikflow-docker-test/sample_1/bam/sample_1.flagstat.txt
/tmp/chikflow-docker-test/sample_1/coverage/sample_1.depth.tsv
/tmp/chikflow-docker-test/sample_1/coverage/sample_1.coverage_summary.csv
/tmp/chikflow-docker-test/sample_1/coverage/sample_1.gene_coverage.csv
/tmp/chikflow-docker-test/sample_1/summary/sample_1.summary.csv
/tmp/chikflow-docker-test/sample_1/assembly/sample_1.consensus.fasta
/tmp/chikflow-docker-test/sample_1/variant_calling/sample_1.variants.vcf.gz
```

## Recommended Next Implementation

1. Add CHIKV genotype/lineage assignment.
   - Define a curated genotype marker/reference strategy.
   - Suggested output: `<sample>/genotyping/*.genotype.csv`.

2. Add phylogeny and final reporting.
   - Build from consensus FASTA outputs and batch summaries.
   - Suggested outputs: phylogeny tree files and `batch_reports/*.html` /
     `batch_reports/*.pdf`.

3. Update docs and tests after each block.
   - Update `README.md`, `docs/output.md`, `docs/architecture.md`, and
     `CHANGELOG.md`.
   - Re-run the focused Singularity test and then Docker once available.

## Useful Commands

Check repository state:

```bash
cd /home/brennovmh/CHIK_FLOW
git status --short --branch
git log --oneline --decorate -5
```

Run the focused Singularity test:

```bash
nextflow run . -profile test,singularity \
  --outdir /tmp/chikflow-singularity-test \
  --skip_fastqc \
  --skip_fastp \
  --skip_multiqc \
  --min_consensus_depth 1
```

Run the Docker test:

```bash
nextflow run . -profile test,docker --outdir /tmp/chikflow-docker-test
```
