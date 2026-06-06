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
  -> select best reference per sample by read k-mer matching
  -> align reads to reference with BWA-MEM
  -> sort/index BAM and collect mapping stats
  -> compute per-base depth and genome coverage summary
  -> summarize coverage across GFF gene/CDS features
  -> call variants and generate masked consensus FASTA
  -> write per-sample mapping/coverage/consensus/variant summary CSV
  -> aggregate per-sample summaries into a batch sample summary CSV
  -> assign nearest-reference genotype/lineage and wild/vaccine source
  -> build batch consensus distance tree with wild/vaccine labels
  -> render batch HTML/PDF report
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
  -> amino-acid variants for CDS-overlapping changes
  -> curated genotype/lineage and wild/vaccine assignment
  -> phylogenetic tree and visualization with wild/vaccine labels
  -> sample and batch PDF/HTML reports with plots
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

The current reporting layer consumes batch summary, genotype CSVs, and Newick
tree output and produces:

- lineage/genotype summary
- HTML report
- PDF report

Future report work should add per-sample coverage plots, batch coverage
heatmaps, richer mutation tables, and rendered tree figures. The implementation
should remain data-driven so reference annotations, lineage panels, and
wild/vaccine source labels can be updated without rewriting report logic.
