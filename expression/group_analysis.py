import pandas as pd
from scipy.stats import ttest_ind
import logging
import os

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def analyze_groups(expression_csv, metadata_csv, output_file):
    logging.info("Loading expression and metadata files...")

    expr = pd.read_csv(expression_csv, index_col=0)
    meta = pd.read_csv(metadata_csv, index_col=0)

    # Align samples
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]

    results = []

    for gene in expr.columns:
        autism = expr.loc[meta["group"] == "autism", gene].dropna()
        control = expr.loc[meta["group"] == "control", gene].dropna()

        if len(autism) < 2 or len(control) < 2:
            logging.warning(f"Not enough samples for {gene}")
            continue

        logfc = autism.mean() - control.mean()
        pval = ttest_ind(autism, control, equal_var=False).pvalue

        results.append({
            "Gene": gene,
            "Autism_Mean": round(autism.mean(), 4),
            "Control_Mean": round(control.mean(), 4),
            "LogFC": round(logfc, 4),
            "P_value": round(pval, 6),
            "Autism_N": len(autism),
            "Control_N": len(control)
        })

    df = pd.DataFrame(results).sort_values("P_value")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, index=False)

    logging.info(f"Group analysis saved to {output_file}")
    return df
