#!/usr/bin/env python
"""Summarize lightweight ATAC-seq QC outputs produced on HPC."""

import argparse
import csv
import re
from pathlib import Path

import yaml


def read_config(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def parse_flagstat(path):
    metrics = {
        "total_reads": None,
        "mapped_reads": None,
        "final_usable_reads": None,
    }
    if not path.exists():
        return metrics

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        count = int(line.split()[0]) if line and line.split()[0].isdigit() else None
        if "in total" in line:
            metrics["total_reads"] = count
        elif " mapped (" in line and "primary" not in line:
            metrics["mapped_reads"] = count
        elif "properly paired" in line:
            metrics["final_usable_reads"] = count
    return metrics


def parse_idxstats(path):
    if not path.exists():
        return None

    total = 0
    mito = 0
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            chrom, _length, mapped, unmapped = line.rstrip("\n").split("\t")[:4]
            mapped = int(mapped)
            unmapped = int(unmapped)
            total += mapped + unmapped
            if chrom in {"chrM", "MT", "M"}:
                mito += mapped + unmapped
    return mito / total if total else None


def parse_markdup(path):
    metrics = {"duplicate_rate": None, "estimated_library_size": None}
    if not path.exists():
        return metrics
    text = path.read_text(encoding="utf-8", errors="ignore")
    read_match = re.search(r"READ:\s*(\d+)", text)
    duplicate_match = re.search(r"DUPLICATE TOTAL:\s*(\d+)", text)
    library_match = re.search(r"ESTIMATED_LIBRARY_SIZE:\s*(\d+)", text)
    if read_match and duplicate_match and int(read_match.group(1)) > 0:
        metrics["duplicate_rate"] = int(duplicate_match.group(1)) / int(read_match.group(1))
    if library_match:
        metrics["estimated_library_size"] = library_match.group(1)
    return metrics


def count_bed_rows(path):
    if not path.exists():
        return None
    count = 0
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def summarize(config_path):
    config = read_config(config_path)
    metadata_path = Path(config["paths"]["metadata"])
    qc_dir = Path(config["paths"]["qc_dir"])
    peak_dir = Path(config["paths"]["peak_dir"])
    summary_path = qc_dir / "atac_qc_summary.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with metadata_path.open(newline="", encoding="utf-8") as handle:
        for sample in csv.DictReader(handle):
            sample_id = sample["sample_id"]
            flagstat = parse_flagstat(qc_dir / "samtools" / f"{sample_id}.filtered.flagstat.txt")
            mito_fraction = parse_idxstats(qc_dir / "samtools" / f"{sample_id}.idxstats.txt")
            markdup = parse_markdup(qc_dir / "samtools" / f"{sample_id}.markdup.metrics.txt")
            peak_count = count_bed_rows(peak_dir / f"{sample_id}_peaks.narrowPeak")
            rows.append(
                {
                    "sample_id": sample_id,
                    "condition": sample.get("condition", ""),
                    "replicate": sample.get("replicate", ""),
                    "total_reads": flagstat["total_reads"],
                    "mapped_reads": flagstat["mapped_reads"],
                    "duplicate_rate": markdup["duplicate_rate"],
                    "mitochondrial_read_fraction": mito_fraction,
                    "final_usable_reads": flagstat["final_usable_reads"],
                    "frip_score": "",
                    "tss_enrichment": "",
                    "fragment_size_distribution": "",
                    "estimated_library_size": markdup["estimated_library_size"],
                    "peak_count": peak_count,
                    "notes": "Blank QC fields require featureCounts/bedtools FRiP, TSS enrichment, or fragment metrics from HPC.",
                }
            )

    fieldnames = [
        "sample_id",
        "condition",
        "replicate",
        "total_reads",
        "mapped_reads",
        "duplicate_rate",
        "mitochondrial_read_fraction",
        "final_usable_reads",
        "frip_score",
        "tss_enrichment",
        "fragment_size_distribution",
        "estimated_library_size",
        "peak_count",
        "notes",
    ]
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"ATAC QC summary saved to {summary_path}")
    return summary_path


def main():
    parser = argparse.ArgumentParser(description="Summarize GSE236992 ATAC-seq QC outputs.")
    parser.add_argument("--config", default="config/atac_config.yaml")
    args = parser.parse_args()
    summarize(args.config)


if __name__ == "__main__":
    main()
