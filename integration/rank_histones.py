import pandas as pd
import numpy as np
import os

EXPR_FILE = "data/processed/autism/autism_vs_control_H1_stats.csv"
DLBC_FILE = "data/processed/dlbc/dlbc_integrated_summary_normalized.csv"

OUTPUT_FILE = "data/processed/integration/ranked_histones.csv"


def normalize(series):
    if series.max() == series.min():
        return series
    return (series - series.min()) / (series.max() - series.min())


def compute_expression_score(df):

    # Replace zero p-values with minimum float
    df["P_value"] = df["P_value"].replace(0, 1e-300)

    # Compute score safely
    df["Expression_raw"] = -np.log10(df["P_value"]) * abs(df["LogFC"])

    # Replace inf and NaN safely
    df["Expression_raw"] = df["Expression_raw"].replace([np.inf, -np.inf], np.nan)
    df["Expression_raw"] = df["Expression_raw"].fillna(0)

    # Normalize
    df["Expression_score"] = normalize(df["Expression_raw"])

    return df



def integrate_scores(expr, dlbc):

    merged = pd.merge(expr, dlbc, on="Gene", how="outer").fillna(0)

    merged = compute_expression_score(merged)

    merged["Mutation_score"] = merged["Mutation_Frequency"]
    merged["CNA_score"] = merged["CNA_Frequency"]

    W_EXPR = 0.4
    W_MUT = 0.3
    W_CNA = 0.3

    merged["Final_score"] = (
        W_EXPR * merged["Expression_score"]
        + W_MUT * merged["Mutation_score"]
        + W_CNA * merged["CNA_score"]
    )

    merged["Final_score"] = merged["Final_score"].replace([np.inf, -np.inf], np.nan)
    merged["Final_score"] = merged["Final_score"].fillna(0)

    merged["Rank"] = (
        merged["Final_score"]
        .rank(ascending=False, method="min")
        .fillna(0)
        .astype(int)
    )

    merged = merged.sort_values("Final_score", ascending=False)

    return merged



def run_ranking():

    expr = pd.read_csv(EXPR_FILE)
    dlbc = pd.read_csv(DLBC_FILE)

    ranked = integrate_scores(expr, dlbc)

    os.makedirs("data/processed/integration", exist_ok=True)

    ranked.to_csv(OUTPUT_FILE, index=False)

    print(f"Ranking saved to {OUTPUT_FILE}")

    return ranked


if __name__ == "__main__":
    run_ranking()
