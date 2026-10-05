# TP53 Variant Landscape: Domain-Level Analysis of ClinVar Missense Variants

A reproducible bioinformatics analysis of how pathogenic and benign missense variants in the tumor suppressor gene **TP53** are distributed across the functional domains of the p53 protein.

## Overview

TP53 is the most frequently mutated gene in human cancer, and germline variants cause Li-Fraumeni syndrome. This project integrates public variant annotations (ClinVar) with protein domain annotations (UniProt) to characterize where clinically significant variants cluster along the protein, and to compare the positional distribution of pathogenic and benign variants.

## Research Question

Are pathogenic TP53 missense variants concentrated in specific functional domains of p53, and how does their positional distribution differ from that of benign variants?

## Data Sources

| Source | Content | Use |
|---|---|---|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Variant classifications and review status | Pathogenic / benign TP53 missense variants |
| [UniProt (P04637)](https://www.uniprot.org/uniprotkb/P04637) | Domain and region annotations | Mapping variants to protein domains |

## Planned Workflow

1. Retrieve TP53 variant records from ClinVar and filter to missense variants with a defined protein change.
2. Parse amino acid positions and classify variants (pathogenic / likely pathogenic vs. benign / likely benign).
3. Map positions to UniProt domain annotations.
4. Quantify and visualize the distribution across domains and along the protein sequence.
5. Summarize findings and discuss limitations (e.g., ClinVar submission bias, review-status heterogeneity).

## Repository Structure

```
tp53-variant-landscape/
├── data/
│   ├── raw/          # Downloaded source files (not tracked)
│   └── processed/    # Cleaned tables
├── notebooks/        # Exploratory and final analysis notebooks
├── scripts/          # Reusable data-processing scripts
├── results/
│   └── figures/      # Output plots
├── requirements.txt
└── README.md
```

## Reproducibility

```bash
git clone https://github.com/Parmida-Mehranfarid/tp53-variant-landscape.git
cd tp53-variant-landscape
pip install -r requirements.txt
```

## Status

In progress. Data retrieval and preprocessing are the first milestone.

## Author

Parmida Mehranfarid, B.Sc. student in Cellular and Molecular Biology
