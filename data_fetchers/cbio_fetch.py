import requests
import pandas as pd
import logging
import os
from tqdm import tqdm
import math

# ---------------------------
# Configuration
# ---------------------------
CBIO_BASE = "https://www.cbioportal.org/api"

# Verified GBM studies with mutation + CNA data
GBM_STUDIES = {
    "TCGA_Nature2008": "gbm_tcga_pub",
    "CPTAC_Cell2021": "gbm_cptac_2021",
    "Columbia_NatMed2019": "gbm_columbia_2019"
}

# Linker histone gene family
H1_GENES = [
    "H1F0",
    "HIST1H1A",
    "HIST1H1B",
    "HIST1H1C",
    "HIST1H1D",
    "HIST1H1E",
    "H1FX"
]

HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

BATCH_SIZE = 150

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# ---------------------------
# API Helpers
# ---------------------------
def get_profiles(study_id):
    """Fetch all molecular profiles for a study"""
    url = f"{CBIO_BASE}/molecular-profiles"
    params = {"studyId": study_id}
    r = requests.get(url, params=params, headers=HEADERS)
    r.raise_for_status()
    return r.json()


def get_mutation_profile(study_id):
    """Resolve mutation molecular profile belonging to study"""
    profiles = get_profiles(study_id)
    for p in profiles:
        pid = p["molecularProfileId"].lower()
        if study_id.lower() in pid and "mutations" in pid:
            return p["molecularProfileId"]
    raise RuntimeError(f"No mutation profile found for {study_id}")


def get_cna_profile(study_id):
    """Resolve CNA/GISTIC molecular profile belonging to study"""
    profiles = get_profiles(study_id)
    for p in profiles:
        pid = p["molecularProfileId"].lower()
        if study_id.lower() in pid and ("cna" in pid or "gistic" in pid):
            return p["molecularProfileId"]
    raise RuntimeError(f"No CNA profile found for {study_id}")


def get_sample_ids(study_id):
    """Fetch sample IDs for a study"""
    url = f"{CBIO_BASE}/studies/{study_id}/samples"
    r = requests.get(url, headers=HEADERS)
    r.raise_for_status()
    return [s["sampleId"] for s in r.json()]


def fetch_mutations_batch(profile_id, genes, sample_batch):
    """POST mutation fetch"""
    url = f"{CBIO_BASE}/molecular-profiles/{profile_id}/mutations/fetch"
    payload = {
        "geneSymbols": genes,
        "sampleIds": sample_batch
    }
    r = requests.post(url, json=payload, headers=HEADERS)
    r.raise_for_status()
    return r.json()


def fetch_cna_batch(profile_id, genes, sample_batch):
    """POST CNA fetch"""
    url = f"{CBIO_BASE}/molecular-profiles/{profile_id}/molecular-data/fetch"
    payload = {
        "geneSymbols": genes,
        "sampleIds": sample_batch
    }
    r = requests.post(url, json=payload, headers=HEADERS)
    r.raise_for_status()
    return r.json()

# ---------------------------
# Main Pipeline
# ---------------------------
def run_gbm_alteration_pipeline(output_dir="data/processed/gbm"):
    os.makedirs(output_dir, exist_ok=True)

    mutation_records = []
    cna_records = []

    for label, study_id in GBM_STUDIES.items():
        logging.info(f"Processing study: {label} ({study_id})")

        try:
            mut_profile = get_mutation_profile(study_id)
            cna_profile = get_cna_profile(study_id)

            logging.info(f"Using mutation profile: {mut_profile}")
            logging.info(f"Using CNA profile: {cna_profile}")

            samples = get_sample_ids(study_id)
            logging.info(f"Samples found: {len(samples)}")

            num_batches = math.ceil(len(samples) / BATCH_SIZE)

            for i in tqdm(range(num_batches), desc=f"{label} batches"):
                batch = samples[i * BATCH_SIZE:(i + 1) * BATCH_SIZE]

                # ---- MUTATIONS ----
                try:
                    mut_data = fetch_mutations_batch(mut_profile, H1_GENES, batch)
                    for r in mut_data:
                        mutation_records.append({
                            "Study": label,
                            "Study_ID": study_id,
                            "Gene": r.get("geneSymbol"),
                            "SampleID": r.get("sampleId"),
                            "Protein_Change": r.get("proteinChange"),
                            "Mutation_Type": r.get("mutationType"),
                            "Variant_Classification": r.get("variantClassification"),
                            "Functional_Impact": r.get("functionalImpactScore"),
                            "Chromosome": r.get("chr"),
                            "Start_Position": r.get("startPosition"),
                            "End_Position": r.get("endPosition")
                        })
                except Exception as e:
                    logging.error(f"Mutation batch {i+1} failed for {label}: {e}")

                # ---- CNA ----
                try:
                    cna_data = fetch_cna_batch(cna_profile, H1_GENES, batch)
                    for r in cna_data:
                        val = r.get("value")
                        if val is None:
                            continue

                        val = int(val)
                        if val != 0:  # Only store altered states
                            cna_records.append({
                                "Study": label,
                                "Study_ID": study_id,
                                "Gene": r.get("geneSymbol"),
                                "SampleID": r.get("sampleId"),
                                "CNA_Value": val
                            })
                except Exception as e:
                    logging.error(f"CNA batch {i+1} failed for {label}: {e}")

        except Exception as e:
            logging.error(f"Failed processing study {label}: {e}")

    # ---------------------------
    # Convert to DataFrames
    # ---------------------------
    mut_df = pd.DataFrame(mutation_records)
    cna_df = pd.DataFrame(cna_records)

    # ---------------------------
    # Save Raw Files
    # ---------------------------
    mut_file = os.path.join(output_dir, "gbm_H1_mutations_raw.csv")
    cna_file = os.path.join(output_dir, "gbm_H1_cna_raw.csv")

    mut_df.to_csv(mut_file, index=False)
    cna_df.to_csv(cna_file, index=False)

    logging.info(f"Saved raw mutation data to {mut_file}")
    logging.info(f"Saved raw CNA data to {cna_file}")

    # ---------------------------
    # Summaries
    # ---------------------------
    mut_summary = (
        mut_df.groupby("Gene")
        .agg(
            Total_Mutations=("SampleID", "count"),
            Unique_Mutated_Samples=("SampleID", "nunique"),
            Studies_With_Mutations=("Study", "nunique"),
            Missense=("Variant_Classification", lambda x: (x == "Missense_Mutation").sum()),
            Truncating=("Variant_Classification", lambda x: x.isin(
                ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins"]
            ).sum())
        .reset_index()
    ) if not mut_df.empty else pd.DataFrame(
        columns=["Gene", "Total_Mutations", "Unique_Mutated_Samples",
                 "Studies_With_Mutations", "Missense", "Truncating"]
    )
    )

    cna_summary = (
        cna_df.groupby("Gene")
        .agg(
            CNA_Events=("SampleID", "count"),
            CNA_Altered_Samples=("SampleID", "nunique"),
            Studies_With_CNA=("Study", "nunique")
        )
        .reset_index()
    ) if not cna_df.empty else pd.DataFrame(
        columns=["Gene", "CNA_Events", "CNA_Altered_Samples", "Studies_With_CNA"]
    )

    # ---------------------------
    # Integrated Summary
    # ---------------------------
    final = pd.merge(mut_summary, cna_summary, on="Gene", how="outer").fillna(0)

    final["Total_Altered_Samples"] = (
        final["Unique_Mutated_Samples"] + final["CNA_Altered_Samples"]
    )

    final = final.sort_values(
        ["Total_Altered_Samples", "Studies_With_CNA", "Studies_With_Mutations"],
        ascending=False
    )

    final_file = os.path.join(output_dir, "gbm_H1_alteration_summary.csv")
    final.to_csv(final_file, index=False)

    logging.info(f"Saved integrated alteration summary to {final_file}")

    return mut_df, cna_df, final


# ---------------------------
# Entry Point
# ---------------------------
if __name__ == "__main__":
    logging.info("Starting GBM H1 alteration pipeline (Mutation + CNA)")
    run_gbm_alteration_pipeline()
