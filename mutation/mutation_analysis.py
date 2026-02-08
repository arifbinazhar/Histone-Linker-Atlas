"""
Mutation Analysis Module

Generates mutation frequency barplot
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------
# PATH CONFIG
# -------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_FILE = PROJECT_ROOT / "data" / "processed" / 'dlbc' /"dlbc_mutation_processed.csv"


PLOTS_DIR = PROJECT_ROOT / "plots"

PLOTS_DIR.mkdir(exist_ok=True)

OUTPUT_FILE = PLOTS_DIR / "dlbc_mutation_barplot.png"


# -------------------------
# LOAD DATA
# -------------------------

def load_mutation():

    df = pd.read_csv(DATA_FILE)

    return df


# -------------------------
# PLOT MUTATION FREQUENCY
# -------------------------

def plot_mutation_frequency(df):

    df_sorted = df.sort_values("Mutation_Count", ascending=False)

    plt.figure(figsize=(8, 5))

    plt.bar(
        df_sorted["Gene"],
        df_sorted["Mutation_Count"]
    )

    plt.title("Mutation Frequency per Histone H1 Gene")

    plt.ylabel("Mutation Count")

    plt.xticks(rotation=45)

    plt.tight_layout()

    plt.savefig(OUTPUT_FILE, dpi=300)

    plt.close()

    print("Saved:", OUTPUT_FILE)


# -------------------------
# MAIN
# -------------------------

def main():

    df = load_mutation()

    plot_mutation_frequency(df)


if __name__ == "__main__":
    main()
