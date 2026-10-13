
"""
Step 1: Clean the ClinVar TP53 export.

Input:
    data/raw/clinvar_tp53_missense.txt

Output:
    data/processed/tp53_missense_clean.csv
"""

from pathlib import Path
import re
import pandas as pd

RAW = Path("data/raw/clinvar_tp53_missense.txt")
OUT = Path("data/processed/tp53_missense_clean.csv")

TRANSCRIPT = "NM_000546.6(TP53)"
PROTEIN_LENGTH = 393

PATHOGENIC = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Pathogenic/Likely pathogenic/Pathogenic, low penetrance",
}

BENIGN = {
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
}

# Three-letter amino acid codes, including Ter for stop.
PROTEIN_CHANGE = re.compile(
    r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)"
)


def classify(label):
    if not isinstance(label, str):
        return "Other"

    label = label.strip()

    if label in PATHOGENIC:
        return "Pathogenic"

    if label in BENIGN:
        return "Benign"

    if label == "Uncertain significance":
        return "VUS"

    if label.startswith("Conflicting"):
        return "Conflicting"

    return "Other"


def main():
    if not RAW.exists():
        raise FileNotFoundError(
            f"Input file not found: {RAW}\n"
            "Download the ClinVar export and place it "
            "at this path."
        )

    df = pd.read_csv(RAW, sep="\t", index_col=False)

    required = {
        "Name",
        "Germline classification",
        "Germline review status",
        "Condition(s)",
        "Accession",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}\n"
            f"Available columns: {list(df.columns)}"
        )

    print(f"Rows in raw file: {len(df)}")

    # Keep records named on the selected TP53 transcript.
    names = df["Name"].fillna("").astype(str)
    df = df[names.str.startswith(TRANSCRIPT)].copy()

    print(f"Rows on selected transcript: {len(df)}")

    # Extract reference amino acid, position and alternate amino acid.
    parsed = df["Name"].str.extract(PROTEIN_CHANGE)
    parsed.columns = ["ref_aa", "position", "alt_aa"]

    df = pd.concat([df, parsed], axis=1)

    no_protein_change = df["position"].isna().sum()
    df = df.dropna(subset=["position"]).copy()

    df["position"] = pd.to_numeric(
        df["position"], errors="coerce"
    )

    invalid_position = (
        df["position"].isna()
        | ~df["position"].between(1, PROTEIN_LENGTH)
        | (df["position"] % 1 != 0)
    )

    print(f"Rows without parsed protein change: {no_protein_change}")
    print(f"Rows with invalid protein position: {invalid_position.sum()}")

    df = df.loc[~invalid_position].copy()
    df["position"] = df["position"].astype(int)

    # Exclude synonymous changes and stop-gain changes.
    # This keeps amino-acid substitutions for the missense analysis.
    is_substitution = (
        df["ref_aa"].ne(df["alt_aa"])
        & df["alt_aa"].ne("Ter")
    )

    excluded_change = (~is_substitution).sum()
    df = df.loc[is_substitution].copy()

    print(f"Excluded non-missense changes: {excluded_change}")
    print(f"Candidate missense records: {len(df)}")

    df["group"] = df["Germline classification"].apply(classify)

    print("\nClassification groups:")
    print(df["group"].value_counts(dropna=False))

    print("\nReview status:")
    print(df["Germline review status"].value_counts(dropna=False))

    # Do not automatically deduplicate records:
    # different ClinVar accessions may contain distinct evidence.
    clean = df[
        [
            "Name",
            "ref_aa",
            "position",
            "alt_aa",
            "Germline classification",
            "Germline review status",
            "group",
            "Condition(s)",
            "Accession",
        ]
    ].rename(
        columns={
            "Germline classification": "classification",
            "Germline review status": "review_status",
            "Condition(s)": "conditions",
        }
    )

    clean = clean.sort_values(
        ["position", "Name"]
    ).reset_index(drop=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    clean.to_csv(OUT, index=False)

    print(f"\nFinal rows: {len(clean)}")
    print(f"Saved to: {OUT}")


if __name__ == "__main__":
    main()
