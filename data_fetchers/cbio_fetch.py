"""
cbio_fetch.py

Fetch CNA and mutation data from cBioPortal GitHub DataHub using sparse checkout,
store raw data in data/raw/,
and processed outputs in data/processed/
"""

import pandas as pd
import subprocess
import shutil
from pathlib import Path
import os
import stat
import time
from data_fetchers.gene_map import H1_GENE_MAP


# -------------------------
# CONFIG
# -------------------------

DATASET = "dlbc_tcga_pan_can_atlas_2018"

REPO_URL = "https://github.com/cBioPortal/datahub.git"

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = PROJECT_ROOT / "data" / "raw" / DATASET
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "dlbc"

TEMP_DIR = PROJECT_ROOT / "datahub_temp"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

H1_GENES = H1_GENE_MAP 


# [
#     "HIST1H1A",
#     "HIST1H1B",
#     "HIST1H1C",
#     "HIST1H1D",
#     "HIST1H1E",
#     "H1F0",
#     "H1FX"
# ]


def safe_rmtree(path, retries=5, delay=1):
    """
    Safely remove directory on Windows by handling file locks.
    """
    def onerror(func, path, exc_info):
        try:
            os.chmod(path, stat.S_IWRITE)
            func(path)
        except Exception:
            pass

    for i in range(retries):
        try:
            shutil.rmtree(path, onerror=onerror)
            return
        except PermissionError:
            time.sleep(delay)

    print(f"Warning: Could not fully remove {path}")


# -------------------------
# STEP 1 — Sparse checkout
# -------------------------

def sparse_clone():

    if RAW_DIR.exists() and any(RAW_DIR.iterdir()):
        print("Dataset already exists:", RAW_DIR)
        return

    print("Cloning required dataset only using sparse checkout...")

    if TEMP_DIR.exists():
        safe_rmtree(TEMP_DIR)

    subprocess.run([
        "git", "clone",
        "--filter=blob:none",
        "--no-checkout",
        REPO_URL,
        str(TEMP_DIR)
    ], check=True)

    subprocess.run(
        ["git", "sparse-checkout", "init", "--cone"],
        cwd=TEMP_DIR,
        check=True
    )

    subprocess.run(
        ["git", "sparse-checkout", "set",
         f"public/{DATASET}"],
        cwd=TEMP_DIR,
        check=True
    )

    subprocess.run(
        ["git", "checkout", "master"],
        cwd=TEMP_DIR,
        check=True
    )

    source = TEMP_DIR / "public" / DATASET

    RAW_DIR.parent.mkdir(parents=True, exist_ok=True)

    shutil.move(str(source), str(RAW_DIR))

    safe_rmtree(TEMP_DIR)

    print("Dataset stored at:", RAW_DIR)



# -------------------------
# STEP 2 — Load data
# -------------------------

def load_data():

    cna_file = RAW_DIR / "data_cna.txt"
    maf_file = RAW_DIR / "data_mutations.txt"

    print("Loading CNA...")
    cna = pd.read_csv(cna_file, sep="\t")

    print("Loading Mutation...")
    maf = pd.read_csv(
        maf_file,
        sep="\t",
        comment="#",
        low_memory=False
    )

    maf["Hugo_Symbol"] = maf["Hugo_Symbol"].astype(str).str.strip()

    return cna, maf


# -------------------------
# STEP 3 — Process CNA
# -------------------------

def process_cna(cna):

    cna = cna.set_index("Hugo_Symbol")

    cna = cna.apply(pd.to_numeric, errors="coerce")

    cna_h1 = cna.loc[cna.index.intersection(H1_GENES)]

    cna_binary = cna_h1.applymap(
        lambda x: 1 if abs(x) >= 1 else 0
    )

    cna_summary = (
        cna_binary.sum(axis=1)
        .reset_index()
    )

    cna_summary.columns = [
        "Gene",
        "CNA_Altered_Samples"
    ]

    return cna_summary


# -------------------------
# STEP 4 — Process mutation
# -------------------------

def process_mutation(maf):

    maf_h1 = maf[maf["Hugo_Symbol"].isin(H1_GENES)]

    mutation_summary = (
        maf_h1.groupby("Hugo_Symbol")
        .agg(
            Mutation_Count=("Tumor_Sample_Barcode", "count"),
            Unique_Patients=("Tumor_Sample_Barcode", "nunique"),
            Missense=("Variant_Classification",
                      lambda x: (x == "Missense_Mutation").sum()),
            Truncating=("Variant_Classification",
                        lambda x: x.isin([
                            "Nonsense_Mutation",
                            "Frame_Shift_Del",
                            "Frame_Shift_Ins",
                            "Splice_Site"
                        ]).sum())
        )
        .reset_index()
        .rename(columns={"Hugo_Symbol": "Gene"})
    )

    return mutation_summary


# -------------------------
# STEP 5 — Integrate
# -------------------------

def integrate(cna_summary, mutation_summary):

    final = pd.merge(
        mutation_summary,
        cna_summary,
        on="Gene",
        how="outer"
    ).fillna(0)

    final["Total_Altered"] = (
        final["Unique_Patients"]
        + final["CNA_Altered_Samples"]
    )

    final = final.sort_values(
        "Total_Altered",
        ascending=False
    )

    return final


# -------------------------
# STEP 6 — Save outputs
# -------------------------

def save_outputs(cna, mutation, integrated):

    cna_file = PROCESSED_DIR / "dlbc_cna_processed.csv"
    mut_file = PROCESSED_DIR / "dlbc_mutation_processed.csv"
    int_file = PROCESSED_DIR / "dlbc_integrated_summary.csv"

    cna.to_csv(cna_file, index=False)
    mutation.to_csv(mut_file, index=False)
    integrated.to_csv(int_file, index=False)

    print("\nSaved files:")
    print(cna_file)
    print(mut_file)
    print(int_file)


# -------------------------
# MAIN
# -------------------------

def main():

    sparse_clone()

    cna, maf = load_data()

    cna_summary = process_cna(cna)

    mutation_summary = process_mutation(maf)

    integrated = integrate(
        cna_summary,
        mutation_summary
    )

    save_outputs(
        cna_summary,
        mutation_summary,
        integrated
    )

    print("\ncbio_fetch pipeline completed successfully.")


if __name__ == "__main__":
    main()
