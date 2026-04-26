#!/usr/bin/env python
"""Link GSE236992 differential ATAC peaks to H1 genes and chromatin regulators."""

from pathlib import Path

import pandas as pd
import yaml


DEFAULT_CONFIG = "config/atac_config.yaml"


def load_config(path=DEFAULT_CONFIG):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _normalize_gene(value):
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def integrate_atac_histones(config_path=DEFAULT_CONFIG):
    config = load_config(config_path)
    da_file = Path(config["paths"]["differential_dir"]) / "differential_accessibility_GSE236992.csv"
    annotation_file = Path(config["paths"]["annotation_dir"]) / "peak_gene_links.csv"
    output_file = Path(config["paths"]["integration_output"])

    if not da_file.exists():
        raise FileNotFoundError(
            f"Missing differential accessibility file: {da_file}. "
            "Run the HPC ATAC workflow before lightweight repo integration."
        )

    da = pd.read_csv(da_file)
    if annotation_file.exists():
        anno = pd.read_csv(annotation_file)
        gene_cols = [col for col in ("peak_id", "nearest_gene", "distance_to_tss", "annotation") if col in anno.columns]
        da = da.merge(anno[gene_cols], on="peak_id", how="left", suffixes=("", "_annotated"))
        if "nearest_gene_annotated" in da.columns:
            da["nearest_gene"] = da["nearest_gene_annotated"].combine_first(da.get("nearest_gene"))
            da = da.drop(columns=["nearest_gene_annotated"])

    histone_genes = {gene.upper() for gene in config["biology"]["linker_histone_genes"]}
    autism_genes = {gene.upper() for gene in config["biology"].get("chromatin_autism_genes", [])}

    da["nearest_gene_normalized"] = da.get("nearest_gene", "").map(_normalize_gene)
    da["is_linker_histone_peak"] = da["nearest_gene_normalized"].isin(histone_genes)
    da["is_chromatin_autism_gene_peak"] = da["nearest_gene_normalized"].isin(autism_genes)

    linked = da[da["is_linker_histone_peak"] | da["is_chromatin_autism_gene_peak"]].copy()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    linked.to_csv(output_file, index=False)

    print(f"ATAC histone/chromatin links saved to {output_file}")
    return linked


if __name__ == "__main__":
    integrate_atac_histones()
