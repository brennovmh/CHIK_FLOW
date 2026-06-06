# Development

## Local Validation

Run lightweight checks before opening a pull request:

```bash
python3 -m py_compile \
  bin/aggregate_sample_summaries.py \
  bin/assign_chikv_genotype.py \
  bin/build_chikv_phylogeny.py \
  bin/calculate_gene_coverage.py \
  bin/render_chikflow_report.py \
  bin/select_reference.py \
  bin/summarize_sample.py \
  bin/validate_genotype_references.py \
  bin/validate_reference_panel.py \
  bin/validate_samplesheet.py \
  bin/vcf_to_aa_mutations.py \
  bin/vcf_to_table.py

python3 tests/test_biological_regression.py

python3 bin/validate_samplesheet.py \
  --input assets/samplesheet.csv \
  --output /tmp/chikflow_validated_samplesheet.csv

python3 bin/validate_reference_panel.py \
  --fasta assets/reference/NC_004162.2.fasta \
  --gff assets/reference/NC_004162.2.gff \
  --output /tmp/chikflow_reference_panel.csv

nextflow run . --help
```

Run the focused Docker smoke test:

```bash
nextflow run . \
  -profile test,docker \
  --outdir /tmp/chikflow-ci-test \
  --skip_fastqc \
  --skip_fastp \
  --skip_multiqc \
  --min_consensus_depth 1
```

This test exercises the core analysis path from bundled FASTQ fixtures through
alignment, coverage, consensus, variants, genotyping, phylogeny, sample
summary, batch summary, and final reports. It intentionally skips FastQC,
fastp, and MultiQC so CI stays focused on CHIKV-specific pipeline behavior.

## CI Scope

The GitHub Actions workflow runs the same lightweight Python checks and focused
Docker smoke test on pushes, pull requests, and manual dispatch. A full
end-to-end run without skips should still be used before tagging releases.
