# CHIKV Genotype Reference Panel

The initial curated genotype panel is stored at:

```text
assets/reference/chikv_genotype_references.fasta
```

Headers use machine-readable metadata:

```text
accession|genotype=<value>|lineage=<value>|source=<wild|vaccine>|country=<value>|year=<value>
```

## Included References

| Accession | Genotype | Lineage label | Source | Notes |
| --- | --- | --- | --- | --- |
| `AY726732.1` | West African | `Senegal_37997` | wild | Senegal 37997 West African reference. |
| `NC_004162.2` | ECSA | `S27_African_prototype` | wild | RefSeq S27 African prototype. |
| `DQ443544.2` | ECSA | `IOL_LR2006_OPY1` | wild | La Reunion 2006 OPY1, Indian Ocean lineage. |
| `KJ451624.1` | Asian | `Caribbean_99659` | wild | Caribbean Asian genotype reference used in Mexico genome report. |
| `KP851709.1` | Asian | `Mexico_InDRE51` | wild | Mexico InDRE51 Asian genotype complete genome. |
| `EF452493.1` | Asian | `AF15561_181_25_parent` | wild | Parent strain for 181/25. |
| `L37661.3` | Asian | `181_25_TSI_GSD_218` | vaccine | Live-attenuated 181/25 vaccine strain. |

## Interpretation

The pipeline classifies each consensus as wild-like, vaccine-like, or unknown
by nearest-reference comparison against this panel. This is a genomic similarity
classification, not a regulatory determination of vaccine product identity.

For stronger vaccine-strain discrimination, add product-specific references
when public complete sequences or validated local surrogate sequences are
available. In particular, VLA1553/IXCHIQ is based on LR2006-OPY1 with an nsP3
deletion, so a validated VLA1553-like reference should be added when available
for operational use.
