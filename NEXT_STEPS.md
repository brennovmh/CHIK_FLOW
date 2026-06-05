# CHIK-FLOW Next Steps

## Current State

The local checkout is synced with GitHub `main`.

Latest implemented blocks:

- Real CHIKV RefSeq reference: `NC_004162.2`
- BWA-MEM alignment
- sorted/indexed BAM and mapping stats
- per-base depth and genome coverage summary
- masked consensus FASTA with `bcftools`
- GFF-derived gene/CDS coverage summary

Last known main commit:

```text
51de63d Add GFF gene coverage summary
```

## Docker Validation Status

Attempted in the resumed session:

```bash
docker run --rm hello-world
```

Result:

```text
permission denied while trying to connect to the docker API at unix:///var/run/docker.sock
```

Open a new terminal/session so the Docker group membership is applied, then run
the same command again.

## After Docker Works

Validate CHIK-FLOW with Docker:

```bash
cd /home/brennovmh/CHIK_FLOW
nextflow run . -profile test,docker --outdir /tmp/chikflow-docker-test
```

If the output directory already exists, use a new one:

```bash
nextflow run . -profile test,docker --outdir /tmp/chikflow-docker-test-2
```

## Implemented After Resume

1. Added a per-sample summary table.
   - Combines mapping stats, genome coverage, GFF feature coverage highlights,
     consensus metrics, low-coverage masking, and VCF record counts.
   - Output: `<sample>/summary/*.summary.csv`.

## Recommended Next Implementation

1. Add variant table conversion.
   - Convert VCF to a readable CSV.
   - Suggested output: `<sample>/variant_calling/*.variants.csv`.

2. Add amino-acid mutation reporting.
   - Use GFF CDS coordinates to translate coding changes.
   - Suggested output: `<sample>/variant_calling/*.aa_mutations.csv`.

3. Add batch-level reporting.
   - Aggregate sample summaries into one CSV.
   - Later use it for HTML/PDF reports.

## Expected Test Outputs

The Docker test run should create files like:

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

## Useful Commands

Check repository state:

```bash
cd /home/brennovmh/CHIK_FLOW
git status --short --branch
git log --oneline --decorate -5
```

Run the focused Singularity test used before Docker validation:

```bash
nextflow run . -profile test,singularity \
  --outdir /tmp/chikflow-singularity-test \
  --skip_fastqc \
  --skip_fastp \
  --skip_multiqc \
  --min_consensus_depth 1
```

Run the full Docker test:

```bash
nextflow run . -profile test,docker --outdir /tmp/chikflow-docker-test
```
