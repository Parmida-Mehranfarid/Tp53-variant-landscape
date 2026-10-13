# 🧬 TP53 Missense Variant Landscape

### A Python-Based Bioinformatics Analysis of Clinical Variant Distribution Across the p53 Protein

[![Python](https://img.shields.io/badge/Python-3.x-blue?logo=python&logoColor=white)](https://www.python.org/)
[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-2E8B57)](https://www.ncbi.nlm.nih.gov/clinvar/)
[![Data Source](https://img.shields.io/badge/Data-NCBI%20ClinVar-orange)](https://www.ncbi.nlm.nih.gov/clinvar/)

## 📌 Project Overview

The **TP53 Missense Variant Landscape** project explores the distribution of missense variants across major functional regions of the human p53 protein and investigates the association between protein location and clinical classifications reported in ClinVar.

The *TP53* gene encodes p53, a tumor suppressor protein involved in cell-cycle regulation, DNA-damage responses, cellular senescence, and apoptosis. Alterations affecting p53 function are frequently associated with human cancer, making TP53 an important subject in cancer biology, molecular genetics, and clinical genomics.

Using Python, this project processes ClinVar-derived variant records, extracts protein-level substitutions, maps variants to predefined protein regions, and applies statistical analysis to compare Pathogenic and Benign classifications.

The primary objective is to demonstrate a reproducible bioinformatics workflow connecting biological data processing, protein-region annotation, data visualization, and statistical interpretation.

## 🎯 Research Questions

- How are TP53 missense protein changes distributed across major regions of the p53 protein?
- What proportion of Pathogenic and Benign protein changes is observed within the DNA-binding domain?
- Is location within the DNA-binding domain associated with Pathogenic rather than Benign classification in the analyzed dataset?

## 🗂️ Data Source

Variant records were obtained from **NCBI ClinVar**, a public archive of reports on the relationships between genetic variants and human health.

**Source:** [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)

The workflow uses protein-level amino-acid substitutions and their corresponding positions to examine the distribution of variants across predefined TP53 protein regions.

The analysis focuses on records assigned to the binary clinical classification groups **Pathogenic** and **Benign**. Other categories, including Variants of Uncertain Significance (VUS), conflicting interpretations, and other classifications, are not included in the binary comparison.

## ⚙️ Methods and Workflow

The project follows a multi-step computational workflow:

1. **Data preprocessing:** Read the downloaded ClinVar-derived dataset and extract relevant variant information.
2. **Protein-change parsing:** Identify amino-acid substitutions and protein positions from the annotated variant descriptions.
3. **Classification filtering:** Retain records classified as Pathogenic or Benign for the primary comparison.
4. **Variant-level quality control:** Identify repeated protein changes and assess whether the same amino-acid change has both Pathogenic and Benign records.
5. **Conflict handling:** Exclude protein changes with conflicting binary classifications from the final comparison.
6. **Protein-region annotation:** Assign each protein change to a predefined TP53 region according to its amino-acid position.
7. **Descriptive analysis:** Calculate Pathogenic and Benign counts and their percentages within the two classification groups for each region.
8. **Statistical testing:** Use Fisher's exact test to assess the association between DNA-binding-domain location and clinical classification.
9. **Visualization and output:** Export a summary table and generate figures illustrating the regional and positional distribution of protein changes.

### Protein regions analyzed

| Protein region | Amino-acid positions |
|---|---:|
| Transactivation | 1–44 |
| N-terminal linker | 45–101 |
| DNA-binding domain | 102–292 |
| Linker | 293–324 |
| Oligomerization | 325–356 |
| C-terminal regulatory region | 357–393 |

These intervals are predefined analytical regions used in this project.

## 📊 Key Results

The final analysis included **388 unique protein changes**, comprising 236 classified as Pathogenic and 152 classified as Benign.

### Classification summary

| Metric | Result |
|---|---:|
| Pathogenic protein changes | 236 |
| Benign protein changes | 152 |
| Total unique protein changes analyzed | 388 |
| Pathogenic changes within the DNA-binding domain | 221 |
| Benign changes within the DNA-binding domain | 42 |
| Pathogenic changes outside the DNA-binding domain | 15 |
| Benign changes outside the DNA-binding domain | 110 |
| Odds ratio | 38.59 |
| Fisher's exact test p-value | 3.3041 × 10⁻⁴⁴ |

### Regional distribution

The DNA-binding domain contained 221 of the 236 Pathogenic protein changes (93.6%) and 42 of the 152 Benign protein changes (27.6%).

The percentages above use the total number of protein changes in each classification group as the denominator. They therefore describe the distribution of Pathogenic or Benign changes across regions, rather than the probability that an individual change in a region is Pathogenic.

### Statistical interpretation

Fisher's exact test compared protein changes located within the DNA-binding domain with those located elsewhere in the protein.

The analysis produced an **odds ratio of 38.59** and a **p-value of 3.3041 × 10⁻⁴⁴**, indicating a strong association between DNA-binding-domain location and Pathogenic classification in this dataset.

This result is consistent with the biological importance of the p53 DNA-binding domain. However, it does not establish that location alone determines pathogenicity, nor does it predict the clinical effect of an individual variant. The findings depend on the selected dataset, clinical classifications, filtering strategy, and predefined region boundaries.

## 📈 Visualizations

### 1. Variant distribution by protein region

![TP53 variants by region](figures/variants_by_region.png)

This figure compares the distribution of Pathogenic and Benign protein changes across the six predefined protein regions.

### 2. Variant distribution along the protein sequence

![TP53 variants along protein](figures/variants_along_protein.png)

This figure visualizes the distribution of analyzed protein changes along the TP53 amino-acid sequence.

## 🧰 Technologies and Tools

- **Python** — scripting and computational analysis
- **Pandas** — data cleaning, filtering, grouping, and tabulation
- **Matplotlib** — data visualization
- **SciPy** — Fisher's exact test
- **NCBI ClinVar** — clinical variant annotation source
- **Git and GitHub** — version control and project documentation

## 📁 Repository Structure

```text
Tp53-variant-landscape/
├── data/
│   ├── raw/
│   │   └── clinvar_tp53_missense.txt
│   └── processed/
│       └── tp53_missense_clean.csv
├── scripts/
│   ├── clean_data.py
│   └── map_domains.py
├── figures/
│   ├── variants_by_region.png
│   └── variants_along_protein.png
├── domain_counts.csv
├── requirements.txt
└── README.md
```

The raw and processed data files are shown to document the intended workflow; their availability in a fresh clone depends on which files have been committed to the repository.

## 🚀 Reproducibility

### Requirements

- Python 3.x
- Packages listed in `requirements.txt`

### Run the analysis

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/Parmida-Mehranfarid/Tp53-variant-landscape.git
cd Tp53-variant-landscape
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the preprocessing and analysis scripts in sequence:

```bash
python scripts/clean_data.py
python scripts/map_domains.py
```

The analysis script generates `domain_counts.csv` and the figures in the `figures/` directory.

**Data requirement:** The raw ClinVar-derived input file must be available at the path expected by the preprocessing script. If it is not included in the repository, obtain the appropriate source data and place it at that location before running the workflow.

## ⚠️ Limitations

- ClinVar classifications are evidence-based annotations that may change as additional evidence becomes available.
- The analysis includes only records assigned to the Pathogenic and Benign categories; other clinical classifications are outside the scope of the primary comparison.
- Repeated records and conflicting binary classifications require careful handling because a protein change may be represented by multiple ClinVar records.
- Protein-region boundaries are predefined for this analysis and may not capture every functional distinction within p53.
- The statistical analysis examines an association in the selected dataset and does not establish causality or provide clinical predictions for individual variants.
- The results may be affected by transcript selection, protein-change parsing, dataset composition, and classification practices.

## 💡 Skills Demonstrated

- Biological data processing and quality control
- Python scripting and Pandas-based analysis
- Protein-level variant parsing
- Protein-region annotation
- Data aggregation and percentage calculations
- Statistical hypothesis testing using Fisher's exact test
- Scientific visualization
- Reproducible research workflows
- Bioinformatics documentation and version control

## 👩‍🔬 Author

**Parmida Mehranfarid**  
B.Sc. Student in Cellular and Molecular Biology  
Research interests: Cancer Biology, Human Genetics, Molecular Genetics, and Bioinformatics

## 📚 Disclaimer

This project was developed for educational and research-training purposes. The results are exploratory and are not intended for clinical diagnosis, genetic counseling, or medical decision-making.
