"""
Expression Analysis Module
Generates:
- Normalized expression matrix
- Heatmap of histone H1 expression
- Barplot per histone gene
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# -------------------------
# PATH CONFIG
# -------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_FILE = PROJECT_ROOT / "data" / "processed" / 'autism' / "autism_vs_control_H1_stats.csv"


PLOTS_DIR = PROJECT_ROOT / "plots"

PLOTS_DIR.mkdir(exist_ok=True)

HEATMAP_FILE = PLOTS_DIR / "histone_expression_heatmap.png"
BARPLOT_FILE = PLOTS_DIR / "histone_expression_barplot.png"


# -------------------------
# LOAD DATA
# -------------------------

def load_expression():

    df = pd.read_csv(DATA_FILE, index_col=0)

    return df


# -------------------------
# NORMALIZATION
# -------------------------

def normalize_expression(df):

    # Z-score normalization
    normalized = (df - df.mean(axis=1).values.reshape(-1, 1)) / df.std(axis=1).values.reshape(-1, 1)

    return normalized


# -------------------------
# HEATMAP
# -------------------------

def plot_heatmap(df):

    plt.figure(figsize=(10, 6))

    sns.heatmap(
        df,
        cmap="coolwarm",
        center=0
    )

    plt.title("Histone H1 Expression Heatmap")

    plt.tight_layout()

    plt.savefig(HEATMAP_FILE, dpi=300)

    plt.close()

    print("Saved:", HEATMAP_FILE)


# -------------------------
# BARPLOT
# -------------------------

def plot_barplot(df):

    mean_expr = df.mean(axis=1)

    plt.figure(figsize=(8, 5))

    mean_expr.sort_values().plot(kind="bar")

    plt.title("Mean Expression per Histone H1 Gene")

    plt.ylabel("Expression")

    plt.tight_layout()

    plt.savefig(BARPLOT_FILE, dpi=300)

    plt.close()

    print("Saved:", BARPLOT_FILE)


# -------------------------
# MAIN
# -------------------------

def main():

    df = load_expression()

    norm = normalize_expression(df)

    plot_heatmap(norm)

    plot_barplot(df)


if __name__ == "__main__":
    main()
