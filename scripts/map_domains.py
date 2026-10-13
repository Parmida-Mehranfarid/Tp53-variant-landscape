
"""
Analyze unique TP53 missense protein changes across protein regions.

Input:
    data/processed/tp53_missense_clean.csv

Outputs:
    domain_counts.csv
    figures/variants_by_region.png
    figures/variants_along_protein.png
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

IN = Path("data/processed/tp53_missense_clean.csv")
OUT_TABLE = Path("domain_counts.csv")
FIG_DIR = Path("figures")

REGIONS = [
    ("Transactivation (1-44)", 1, 44),
    ("N-terminal linker (45-101)", 45, 101),
    ("DNA-binding domain (102-292)", 102, 292),
    ("Linker (293-324)", 293, 324),
    ("Oligomerization (325-356)", 325, 356),
    ("C-terminal regulatory (357-393)", 357, 393),
]

COLORS = {"Pathogenic": "#c0392b", "Benign": "#2e86c1"}


def assign_region(position):
    for name, start, end in REGIONS:
        if start <= position <= end:
            return name
    return None


def main():
    if not IN.exists():
        raise FileNotFoundError(
            f"Processed input not found: {IN}\n"
            "Run scripts/clean_data.py first."
        )

    df = pd.read_csv(IN)

    required = {"ref_aa", "position", "alt_aa", "group"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = df[df["group"].isin(["Pathogenic", "Benign"])].copy()

    df["variant_key"] = (
        df["ref_aa"].astype(str)
        + df["position"].astype(str)
        + df["alt_aa"].astype(str)
    )

    # Keep a protein change only if all its retained records agree
    # on one classification. Exclude classification conflicts.
    classifications_per_variant = df.groupby("variant_key")["group"].nunique()
    conflicting_keys = classifications_per_variant[
        classifications_per_variant > 1
    ].index

    print(f"Pathogenic/benign records before deduplication: {len(df)}")
    print(f"Unique protein changes before conflict filtering: "
          f"{df['variant_key'].nunique()}")
    print(f"Protein changes with conflicting binary classifications: "
          f"{len(conflicting_keys)}")

    df = df[~df["variant_key"].isin(conflicting_keys)].copy()

    # At this point each protein change has one classification.
    df = df.drop_duplicates(subset="variant_key", keep="first").copy()

    print(f"Unique protein changes included: {len(df)}")
    print("\nClassification counts:")
    print(df["group"].value_counts().to_string())

    df["region"] = df["position"].apply(assign_region)

    if df["region"].isna().any():
        unmapped = df.loc[df["region"].isna(), ["variant_key", "position"]]
        print("\nUnmapped protein changes:")
        print(unmapped.to_string(index=False))
        raise ValueError("Some protein changes could not be mapped.")

    order = [region[0] for region in REGIONS]

    counts = (
        df.groupby(["region", "group"])
        .size()
        .unstack(fill_value=0)
        .reindex(order, fill_value=0)
    )

    for group in ["Pathogenic", "Benign"]:
        if group not in counts.columns:
            counts[group] = 0

    counts = counts[["Pathogenic", "Benign"]]
    pct = counts.div(counts.sum(axis=0).replace(0, float("nan")), axis=1) * 100
    pct.columns = ["Pathogenic_%", "Benign_%"]

    table = pd.concat([counts, pct.round(1)], axis=1)
    table.to_csv(OUT_TABLE, index_label="region")

    print("\nCounts and percentages by region:")
    print(table.to_string())

    dbd = "DNA-binding domain (102-292)"
    p_in = int(counts.loc[dbd, "Pathogenic"])
    b_in = int(counts.loc[dbd, "Benign"])
    p_out = int(counts["Pathogenic"].sum() - p_in)
    b_out = int(counts["Benign"].sum() - b_in)

    print("\nDNA-binding domain vs rest of protein")
    print("Fisher's exact test")

    if min(p_in, b_in, p_out, b_out) == 0:
        print("A zero cell is present; interpret the odds ratio carefully.")

    odds_ratio, p_value = fisher_exact(
        [[p_in, b_in], [p_out, b_out]],
        alternative="two-sided",
    )

    print(f"Pathogenic: {p_in} in domain / {p_out} outside")
    print(f"Benign:     {b_in} in domain / {b_out} outside")
    print(f"Odds ratio: {odds_ratio:.4g}")
    print(f"P-value:    {p_value:.4e}")

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    x = list(range(len(order)))
    width = 0.38

    ax.bar(
        [i - width / 2 for i in x],
        pct["Pathogenic_%"].fillna(0),
        width,
        label="Pathogenic",
        color=COLORS["Pathogenic"],
    )
    ax.bar(
        [i + width / 2 for i in x],
        pct["Benign_%"].fillna(0),
        width,
        label="Benign",
        color=COLORS["Benign"],
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [name.replace(" (", "\n(") for name in order],
        fontsize=8,
    )
    ax.set_ylabel("Share of unique variants in group (%)")
    ax.set_title("Unique TP53 missense variants by protein region")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "variants_by_region.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4))
    shades = ["#f4f6f7", "#e5e8e8"]

    for i, (_, start, end) in enumerate(REGIONS):
        ax.axvspan(
            start - 0.5,
            end + 0.5,
            color=shades[i % 2],
            zorder=0,
        )

    for group in ["Pathogenic", "Benign"]:
        per_position = (
            df[df["group"] == group]
            .groupby("position")["variant_key"]
            .nunique()
        )
        values = (
            per_position.values
            if group == "Pathogenic"
            else -per_position.values
        )
        ax.bar(
            per_position.index,
            values,
            width=1.0,
            color=COLORS[group],
            label=group,
            zorder=2,
        )

    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_xlim(0, 394)
    ax.set_xlabel("Amino acid position")
    ax.set_ylabel("Unique variants per position\n(Benign shown below axis)")
    ax.set_title("Unique TP53 missense variants along the protein")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "variants_along_protein.png", dpi=200)
    plt.close(fig)

    print(f"\nTable saved to: {OUT_TABLE}")
    print(f"Figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()
