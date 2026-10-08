"""
Step 2: Map TP53 missense variants to p53 protein regions and compare
pathogenic vs. benign distributions.

Input : data/processed/tp53_missense_clean.csv   (output of clean_data.py)
Output: results/domain_counts.csv
        results/figures/variants_by_region.png
        results/figures/variants_along_protein.png

Region boundaries (UniProt P04637, human p53, 393 aa)
-----------------------------------------------------
- 1-44    Transcription activation (acidic)   [UniProt 'Region']
- 102-292 DNA-binding domain                  [UniProt 'DNA binding']
- 325-356 Oligomerization (tetramerization)   [UniProt 'Region']
The stretches between these annotated regions (45-101, 293-324, 357-393)
are grouped here as linker / flanking regions for the purpose of this
analysis. They are an analysis choice, not UniProt domains.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

IN = Path("data/processed/tp53_missense_clean.csv")
OUT_TABLE = Path("results/domain_counts.csv")
FIG_DIR = Path("results/figures")

REGIONS = [
    ("Transactivation (1-44)", 1, 44),
    ("N-terminal linker (45-101)", 45, 101),
    ("DNA-binding domain (102-292)", 102, 292),
    ("Linker (293-324)", 293, 324),
    ("Oligomerization (325-356)", 325, 356),
    ("C-terminal regulatory (357-393)", 357, 393),
]
COLORS = {"Pathogenic": "#c0392b", "Benign": "#2e86c1"}


def assign_region(pos):
    for name, start, end in REGIONS:
        if start <= pos <= end:
            return name
    return None


def main():
    df = pd.read_csv(IN)
    df = df[df["group"].isin(["Pathogenic", "Benign"])].copy()
    df["region"] = df["position"].apply(assign_region)
    print(f"Variants analysed: {len(df)}")
    print(df["group"].value_counts().to_string(), "\n")

    order = [r[0] for r in REGIONS]
    counts = (
        df.groupby(["region", "group"]).size().unstack(fill_value=0).reindex(order, fill_value=0)
    )
    for g in ["Pathogenic", "Benign"]:
        if g not in counts.columns:
            counts[g] = 0
    counts = counts[["Pathogenic", "Benign"]]
    pct = counts.div(counts.sum(axis=0), axis=1) * 100
    pct.columns = ["Pathogenic_%", "Benign_%"]
    table = pd.concat([counts, pct.round(1)], axis=1)
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_TABLE)
    print(table.to_string(), "\n")

    # Fisher's exact test: DNA-binding domain vs. rest of protein
    dbd = "DNA-binding domain (102-292)"
    p_in = counts.loc[dbd, "Pathogenic"]
    b_in = counts.loc[dbd, "Benign"]
    p_out = counts["Pathogenic"].sum() - p_in
    b_out = counts["Benign"].sum() - b_in
    odds, pval = fisher_exact([[p_in, b_in], [p_out, b_out]])
    print("DNA-binding domain vs rest of protein (Fisher's exact test)")
    print(f"  Pathogenic: {p_in} in DBD / {p_out} outside")
    print(f"  Benign:     {b_in} in DBD / {b_out} outside")
    print(f"  Odds ratio = {odds:.2f}, p-value = {pval:.2e}\n")

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Figure 1: percentage of each group per region
    fig, ax = plt.subplots(figsize=(9, 5))
    x = range(len(order))
    w = 0.38
    ax.bar([i - w / 2 for i in x], pct["Pathogenic_%"], w, label="Pathogenic", color=COLORS["Pathogenic"])
    ax.bar([i + w / 2 for i in x], pct["Benign_%"], w, label="Benign", color=COLORS["Benign"])
    ax.set_xticks(list(x))
    ax.set_xticklabels([o.replace(" (", "\n(") for o in order], fontsize=8)
    ax.set_ylabel("Share of variants in group (%)")
    ax.set_title("TP53 missense variants by protein region (ClinVar)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(FIG_DIR / "variants_by_region.png", dpi=200)
    plt.close()

    # Figure 2: variants along the protein sequence
    fig, ax = plt.subplots(figsize=(11, 4))
    shade = ["#f4f6f7", "#e5e8e8"]
    for i, (name, start, end) in enumerate(REGIONS):
        ax.axvspan(start - 0.5, end + 0.5, color=shade[i % 2], zorder=0)
    for g in ["Pathogenic", "Benign"]:
        per_pos = df[df["group"] == g].groupby("position").size()
        ax.bar(per_pos.index, per_pos.values if g == "Pathogenic" else -per_pos.values,
               width=1.0, color=COLORS[g], label=g, zorder=2)
    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_xlim(0, 394)
    ax.set_xlabel("Amino acid position")
    ax.set_ylabel("Variants per position\n(benign shown below axis)")
    ax.set_title("Distribution of TP53 missense variants along the p53 protein")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "variants_along_protein.png", dpi=200)
    plt.close()
    print(f"Figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
