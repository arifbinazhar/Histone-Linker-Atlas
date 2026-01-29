import os
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def extract_groups_from_gse(gse, output_dir, disease_name="autism"):
    """
    Extract sample phenotype labels from GEO metadata.
    Attempts to classify samples into disease vs control based on sample characteristics.
    """
    records = []

    for gsm_name, gsm in gse.gsms.items():
        group = "unknown"

        # Search in sample metadata for disease/control keywords
        meta_text = " ".join(
            [" ".join(map(str, v)) for v in gsm.metadata.values()]
        ).lower()

        if disease_name.lower() in meta_text:
            group = "autism"
        elif "control" in meta_text or "healthy" in meta_text or "normal" in meta_text:
            group = "control"

        records.append({
            "SampleID": gsm_name,
            "group": group
        })

    df = pd.DataFrame(records).set_index("SampleID")

    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, "metadata_autism_groups.csv")
    df.to_csv(out_file)

    logging.info(f"Saved sample metadata to {out_file}")
    return out_file
