import pandas as pd
import os

INPUT_FILE = "data/processed/dlbc/dlbc_integrated_summary.csv"
OUTPUT_FILE = "data/processed/dlbc/dlbc_integrated_summary_normalized.csv"

# Replace this with actual cohort size if known
# DLBC cohort is typically 48 patients in TCGA DLBC
TOTAL_PATIENTS = 48


def normalize_dlbc():

    df = pd.read_csv(INPUT_FILE)

    # Normalize mutation recurrence
    df["Mutation_Frequency"] = df["Unique_Patients"] / TOTAL_PATIENTS

    # Normalize CNA recurrence
    df["CNA_Frequency"] = df["CNA_Altered_Samples"] / TOTAL_PATIENTS

    # Normalize overall alteration
    df["Alteration_Frequency"] = df["Total_Altered"] / TOTAL_PATIENTS

    # Optional: normalize mutation impact score
    df["Functional_Mutation_Frequency"] = (
        df["Missense"] + df["Truncating"]
    ) / TOTAL_PATIENTS

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"Normalized DLBC summary saved to {OUTPUT_FILE}")

    return df


if __name__ == "__main__":
    normalize_dlbc()
