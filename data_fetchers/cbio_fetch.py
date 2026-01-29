import requests
import pandas as pd
import logging
import os
from tqdm import tqdm
import math

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

CBIO_BASE = "https://www.cbioportal.org/api"

# Explicit TCGA Glioblastoma study
GBM_STUDY = "gbm_tcga"

H1_GENES = [
    "H1F0",
    "HIST1H1A",
    "HIST1H1B",
    "HIST1H1C",
    "HIST1H1D",
    "HIST1H1E",
    "H1FX"
]

HEADERS = {"Content-Type": "application/json"}


# ---------------------------
# API Helpers
# ---------------------------
def get_mutation_profile(study_id):
    url = f"{CBIO_BASE}/molecular-profiles"
    params = {"studyId": study_id}
    r = requests.get(url, params=params)
    r.raise_for_status()

    profiles = r.json()
    for p in profiles:
        pid = p["molecularProfileId"].lower()
        if "mutations" in pid and study_id in pid:
            return p["molecularProfileId"]

    raise RuntimeError("No GBM mutation profile found.")


def get_sample_ids(study_id):
    url = f"{CBIO_BASE}/studies/{study_id}/samples"
    r = requests.get(url)
    r.raise_for_status()
    samples = r.json()
    return [s["sampleId"] for s in samples]


def fetch_mutations_batch(profile_id, gene, sample_batch):
    url = f"{CBIO_BASE}/molecular-profiles/{profile_id}/mutations/fetch"
    payload = {
        "geneSymbols": [gene],
        "sampleIds": sample_batch
    }

    r = requests.post(url, json=payload, headers=HEADERS)
    r.raise_for_status()
    return r.json()


# ---------------------------
# Main Pipeline
# ---------------------------
def run_gbm_mutation_pipeline(output_dir="data/processed/gbm"):
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Resolving GBM mutation profile...")
    profile_id = get_mutation_profile(GBM_STUDY)
    logging.info(f"Using mutation profile: {profile_id}")

    logging.info("Fetching GBM sample IDs...")
    samples = get_sample_ids(GBM_STUDY)
    logging.info(f"Total samples: {len(samples)}")

    all_records = []
    batch_size = 200
    num_batches = math.ceil(len(samples) / batch_size)

    for gene in tqdm(H1_GENES, desc="Fetching H1 mutations"):
        try:
            for i in range(num_batches):
                batch = samples[i * batch_size:(i + 1) * batch_size]
                records = fetch_mutations_batch(profile_id, gene, batch)

                for r in records:
                    all_records.append({
                        "Gene": gene,
                        "SampleID": r.get("sampleId"),
                        "Protein_Change": r.get("proteinChange"),
                        "Mutation_Type": r.get("mutationType"),
                        "Functional_Impact": r.get("functionalImpactScore")
                    })

        except Exception as e:
            logging.error(f"Failed fetching {gene}: {e}")

    df = pd.DataFrame(all_records)

    if df.empty:
        logging.warning("No mutation records found for H1 genes in GBM.")
        return None, None

    raw_file = os.path.join(output_dir, "gbm_H1_mutations_raw.csv")
    df.to_csv(raw_file, index=False)
    logging.info(f"Saved raw mutation data to {raw_file}")

    summary = (
        df.groupby("Gene")
        .agg(
            Mutation_Count=("SampleID", "count"),
            Unique_Samples=("SampleID", "nunique")
        )
        .reset_index()
        .sort_values("Mutation_Count", ascending=False)
    )

    summary_file = os.path.join(output_dir, "gbm_H1_mutation_summary.csv")
    summary.to_csv(summary_file, index=False)
    logging.info(f"Saved mutation summary to {summary_file}")

    return df, summary


if __name__ == "__main__":
    run_gbm_mutation_pipeline()
