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
        fq1_value = row.get("fastq_1", "").strip()
        fq2_value = row.get("fastq_2", "").strip()
        if not fq1_value or not fq2_value:
            print(f"Skipping {sample}: paired FASTQ paths are incomplete.")
            continue
        fq1 = Path(fq1_value)
        fq2 = Path(fq2_value)
        if not fq1.exists() or not fq2.exists():
            raise FileNotFoundError(
                f"{sample}: missing input FASTQ(s): {fq1} {fq2}. "
                "Check metadata/atac_samples_GSE236992.csv and the download directory."
            )
        out1 = outdir / f"{sample}.trimmed.R1.fastq.gz"
        out2 = outdir / f"{sample}.trimmed.R2.fastq.gz"
        html = outdir / f"{sample}.fastp.html"
        json = outdir / f"{sample}.fastp.json"
        if out1.exists() and out2.exists() and out1.stat().st_size > 0 and out2.stat().st_size > 0:
            print(f"Skipping {sample}: trimmed FASTQs already exist.", flush=True)
            continue
        print(f"Trimming {sample}: {fq1} {fq2}", flush=True)
        subprocess.run(
            [
                "fastp",
                "-i",
                str(fq1),
                "-I",
                str(fq2),
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
        print(f"Finished {sample}: {out1} {out2}", flush=True)
PY
