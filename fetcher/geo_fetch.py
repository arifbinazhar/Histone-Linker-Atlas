import GEOparse
import pandas as pd
import os
import logging
from tqdm import tqdm
from fetcher.gene_map import H1_GENE_MAP

# -------------------------
# Logging Setup
# -------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# -------------------------
# Utilities
# -------------------------
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

# -------------------------
# GEO Dataset Search (Manual List for Stability)
# -------------------------
AUTISM_GEO_DATASETS = [
    "GSE28521",  # Brain tissue autism vs control
    "GSE42133"   # iPSC-derived neurons autism
]

# -------------------------
# Expression Extraction
# -------------------------
def extract_h1_expression(gse, output_dir):
    logging.info(f"Processing {gse.name}")

    matrix = gse.pivot_samples('VALUE')
    matrix.index = matrix.index.str.upper()

    matched_genes = {}
    for official, aliases in H1_GENE_MAP.items():
        for alias in aliases:
            alias = alias.upper()
            if alias in matrix.index:
                matched_genes[official] = matrix.loc[alias]
                break

    if not matched_genes:
        logging.warning(f"No H1 genes found in {gse.name}")
        return None

    df = pd.DataFrame(matched_genes)
    df.index.name = "SampleID"

    out_file = os.path.join(output_dir, f"{gse.name}_H1_expression.csv")
    df.to_csv(out_file)

    logging.info(f"Saved H1 expression to {out_file}")
    return df

# -------------------------
# Main Pipeline
# -------------------------
def run_autism_pipeline(output_dir="data/processed/autism"):
    ensure_dir(output_dir)

    all_results = []

    for gse_id in tqdm(AUTISM_GEO_DATASETS, desc="Downloading GEO datasets"):
        try:
            gse = GEOparse.get_GEO(gse_id, destdir="data/raw")
            df = extract_h1_expression(gse, output_dir)

            if df is not None:
                df["Dataset"] = gse_id
                all_results.append(df)

        except Exception as e:
            logging.error(f"Failed processing {gse_id}: {e}")

    if all_results:
        merged = pd.concat(all_results)
        merged_file = os.path.join(output_dir, "merged_autism_H1_expression.csv")
        merged.to_csv(merged_file)
        logging.info(f"Merged dataset saved to {merged_file}")

    logging.info("Autism GEO pipeline completed.")

if __name__ == "__main__":
    run_autism_pipeline()
