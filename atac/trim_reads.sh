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
outdir = Path(config["paths"]["trimmed_dir"])
threads = str(config.get("tools", {}).get("threads", 8))
outdir.mkdir(parents=True, exist_ok=True)

with metadata.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        sample = row["sample_id"]
        fq1 = row.get("fastq_1", "").strip()
        fq2 = row.get("fastq_2", "").strip()
        if not fq1 or not fq2:
            print(f"Skipping {sample}: paired FASTQ paths are incomplete.")
            continue
        out1 = outdir / f"{sample}.trimmed.R1.fastq.gz"
        out2 = outdir / f"{sample}.trimmed.R2.fastq.gz"
        html = outdir / f"{sample}.fastp.html"
        json = outdir / f"{sample}.fastp.json"
        subprocess.run(
            [
                "fastp",
                "-i",
                fq1,
                "-I",
                fq2,
                "-o",
                str(out1),
                "-O",
                str(out2),
                "--thread",
                threads,
                "--detect_adapter_for_pe",
                "--html",
                str(html),
                "--json",
                str(json),
            ],
            check=True,
        )
PY
