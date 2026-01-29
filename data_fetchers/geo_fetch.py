import GEOparse
import pandas as pd
import os
import logging
from tqdm import tqdm
from data_fetchers.gene_map import H1_GENE_MAP
from data_fetchers.geo_metadata import extract_groups_from_gse



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

    # Get expression matrix (rows = probes, columns = samples)
    matrix = gse.pivot_samples('VALUE')

    # Get platform annotation (probe → gene symbol)
    gpl = list(gse.gpls.values())[0]
    gpl_table = gpl.table

    # Normalize column names
    gpl_table.columns = [c.upper() for c in gpl_table.columns]

    # Find gene symbol column
    gene_col = None
    for col in gpl_table.columns:
        if "GENE" in col or "SYMBOL" in col:
            gene_col = col
            break

    if gene_col is None:
        logging.warning(f"No gene symbol column found in platform for {gse.name}")
        return None

    # Build probe → gene mapping
    probe_to_gene = dict(zip(gpl_table["ID"], gpl_table[gene_col]))

    # Map probes to genes
    matrix["GENE_SYMBOL"] = matrix.index.map(lambda x: probe_to_gene.get(x, "").upper())

    # Reverse mapping: gene → expression (mean if multiple probes)
    gene_matrix = matrix.groupby("GENE_SYMBOL").mean(numeric_only=True)

    matched_genes = {}
    for official, aliases in H1_GENE_MAP.items():
        for alias in aliases:
            alias = alias.upper()
            if alias in gene_matrix.index:
                matched_genes[official] = gene_matrix.loc[alias]
                break

    if not matched_genes:
        logging.warning(f"No H1 genes found in {gse.name} after probe mapping")
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
            gse = GEOparse.get_GEO(gse_id, destdir="data/raw", how="full", annotate_gpl=True)
            df = extract_h1_expression(gse, output_dir)

            if df is not None:
                extract_groups_from_gse(gse, output_dir)
                df["Dataset"] = gse_id
                all_results.append(df)

            
            # df = extract_h1_expression(gse, output_dir)

            # if df is not None:
            #     df["Dataset"] = gse_id
            #     all_results.append(df)

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
