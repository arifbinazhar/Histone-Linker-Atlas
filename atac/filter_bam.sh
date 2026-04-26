#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config/atac_config.yaml}"

python - "$CONFIG" <<'PY'
import csv
import subprocess
import sys
from pathlib import Path

import yaml

config = yaml.safe_load(open(sys.argv[1], "r", encoding="utf-8"))
metadata = Path(config["paths"]["metadata"])
aligned_dir = Path(config["paths"]["aligned_dir"])
filtered_dir = Path(config["paths"].get("filtered_dir", "results/atac/GSE236992/filtered"))
qc_dir = Path(config["paths"]["qc_dir"])
blacklist = config["reference"]["blacklist_bed"]
mapq = str(config.get("tools", {}).get("mapq", 30))
threads = str(config.get("tools", {}).get("threads", 8))
filtered_dir.mkdir(parents=True, exist_ok=True)
(qc_dir / "samtools").mkdir(parents=True, exist_ok=True)

with metadata.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        sample = row["sample_id"]
        bam = aligned_dir / f"{sample}.sorted.bam"
        primary = filtered_dir / f"{sample}.primary.bam"
        no_mito = filtered_dir / f"{sample}.no_mito.bam"
        name_sorted = filtered_dir / f"{sample}.namesort.bam"
        fixmate = filtered_dir / f"{sample}.fixmate.bam"
        position_sorted = filtered_dir / f"{sample}.positionsort.bam"
        dedup = filtered_dir / f"{sample}.dedup.bam"
        final = filtered_dir / f"{sample}.filtered.bam"
        metrics = qc_dir / "samtools" / f"{sample}.markdup.metrics.txt"

        subprocess.run(["samtools", "view", "-@", threads, "-b", "-F", "1804", "-q", mapq, str(bam), "-o", str(primary)], check=True)
        subprocess.run(["samtools", "idxstats", str(primary)], check=True, stdout=(qc_dir / "samtools" / f"{sample}.idxstats.txt").open("w"))
        cmd_no_mito = f"samtools view -@ {threads} -h {primary} | awk '$3 != \"chrM\" && $3 != \"MT\" || $1 ~ /^@/' | samtools sort -@ {threads} -o {no_mito} -"
        subprocess.run(cmd_no_mito, shell=True, check=True)
        subprocess.run(["samtools", "sort", "-n", "-@", threads, "-o", str(name_sorted), str(no_mito)], check=True)
        subprocess.run(["samtools", "fixmate", "-m", str(name_sorted), str(fixmate)], check=True)
        subprocess.run(["samtools", "sort", "-@", threads, "-o", str(position_sorted), str(fixmate)], check=True)
        subprocess.run(["samtools", "markdup", "-r", "-s", str(position_sorted), str(dedup)], check=True, stderr=metrics.open("w"))
        subprocess.run(["bedtools", "intersect", "-v", "-abam", str(dedup), "-b", blacklist], check=True, stdout=final.open("wb"))
        subprocess.run(["samtools", "index", str(final)], check=True)
        subprocess.run(["samtools", "flagstat", str(final)], check=True, stdout=(qc_dir / "samtools" / f"{sample}.filtered.flagstat.txt").open("w"))
PY
