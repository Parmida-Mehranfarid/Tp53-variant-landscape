"""
Step 1: Clean the ClinVar TP53 export.

Input : data/raw/clinvar_tp53_missense.txt  (ClinVar tabular download)
Output: data/processed/tp53_missense_clean.csv

What this script does
---------------------
1. Reads the ClinVar table.
2. Extracts the amino acid position from the 'Name' column, which uses the
   canonical TP53 transcript (NM_000546.6, 393 aa). The 'Protein change'
   column is NOT used because it lists several isoform numberings at once.
3. Keeps only variants with a defined missense protein change.
4. Groups germline classifications into 'Pathogenic' and 'Benign'.
   Uncertain significance and conflicting classifications are kept but
   labelled separately, and are excluded from the main comparison.
"""

from pathlib import Path
import re
import pandas as pd

RAW = Path("data/raw/clinvar_tp53_missense.txt")
OUT = Path("data/processed/tp53_missense_clean.csv")

PATHOGENIC = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Pathogenic/Likely pathogenic/Pathogenic, low penetrance",
}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}

# Matches e.g. "(p.Arg273Cys)" and captures reference aa, position, alternate aa
PROTEIN_CHANGE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")


def classify(label):
    if label in PATHOGENIC:
        return "Pathogenic"
    if label in BENIGN:
        return "Benign"
    if label == "Uncertain significance":
        return "VUS"
    if isinstance(label, str) and label.startswith("Conflicting"):
        return "Conflicting"
    return "Other"


def main():
    # index_col=False: data rows end with an extra tab, which would otherwise
    # shift all columns by one.
    df = pd.read_csv(RAW, sep="\t", index_col=False)
    print(f"Rows in raw file: {len(df)}")

    # Keep only the canonical TP53 transcript. The export also contains a few
    # variants named on overlapping genes (e.g. WRAP53), whose protein
    # positions are not TP53 positions.
    df = df[df["Name"].str.startswith("NM_000546.6(TP53)")].copy()
    print(f"Rows on canonical TP53 transcript (NM_000546.6): {len(df)}")

    parsed = df["Name"].str.extract(PROTEIN_CHANGE)
    parsed.columns = ["ref_aa", "position", "alt_aa"]
    df = pd.concat([df, parsed], axis=1)

    # Drop intronic/UTR rows without a protein change, and synonymous/nonsense edge cases
    df = df.dropna(subset=["position"]).copy()
    df = df[(df["alt_aa"] != df["ref_aa"]) & (df["alt_aa"] != "Ter")]
    df["position"] = df["position"].astype(int)
    print(f"Missense variants with a protein position: {len(df)}")

    df["group"] = df["Germline classification"].apply(classify)

    clean = df[
        [
            "Name", "ref_aa", "position", "alt_aa",
            "Germline classification", "Germline review status",
            "group", "Condition(s)", "Accession",
        ]
    ].rename(
        columns={
            "Germline classification": "classification",
            "Germline review status": "review_status",
            "Condition(s)": "conditions",
        }
    )
    clean = clean.sort_values("position").reset_index(drop=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    clean.to_csv(OUT, index=False)

    print("\nVariants per group:")
    print(clean["group"].value_counts())
    print(f"\nSaved to {OUT}")


if __name__ == "__main__":
    main()
