# CHIK-FLOW Architecture

CHIK-FLOW follows a Nextflow DSL2 structure inspired by nf-core and RSVrecon.
The first milestone is a clean executable foundation. CHIKV-specific
interpretation will be added as separate modules with test data.

## Current Workflow

```text
samplesheet
  -> validate samplesheet
  -> group and merge lanes per sample
  -> FastQC pre-trim
  -> fastp trimming/filtering
  -> FastQC post-trim
  -> MultiQC
  -> validate reference panel
  -> align reads to reference with BWA-MEM
  -> sort/index BAM and collect mapping stats
  -> compute per-base depth and genome coverage summary
  -> summarize coverage across GFF gene/CDS features
  -> call variants and generate masked consensus FASTA
  -> write per-sample mapping/coverage/consensus/variant summary CSV
```

## Target Workflow

```text
FASTQ
  -> QC and trimming
  -> CHIKV reference-panel screening
  -> best reference selection
  -> alignment
  -> BAM sort/index/stats
  -> per-base depth and genome coverage
  -> per-gene coverage
  -> consensus FASTA
  -> nucleotide variants and readable variant CSV table
  -> amino-acid variants by ORF/protein
  -> genotype/lineage assignment
  -> phylogenetic tree
  -> sample and batch PDF/HTML reports
```

## CHIKV Regions for Coverage Reporting

The planned coverage summary should report whole-genome coverage plus gene or
region coverage for:

- 5'UTR
- nsP1
- nsP2
- nsP3
- nsP4
- C
- E3
- E2
- 6K/TF
- E1
- 3'UTR

Coordinates will come from the selected reference GFF rather than hard-coded
global positions.

## Report Strategy

The reporting layer will consume a per-sample manifest and produce:

- `Report.csv`
- `Consensus.fasta`
- per-sample coverage plots
- batch coverage heatmap
- mutation tables
- lineage/genotype summary
- phylogenetic tree figures
- HTML report
- PDF report

The implementation should remain data-driven so reference annotations and
lineage panels can be updated without rewriting report logic.
