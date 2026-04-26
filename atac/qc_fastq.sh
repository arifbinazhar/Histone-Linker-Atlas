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
fastqc_dir = Path(config["paths"]["qc_dir"]) / "fastqc"
multiqc_dir = Path(config["paths"]["qc_dir"]) / "multiqc"
threads = str(config.get("tools", {}).get("threads", 8))
fastqc_dir.mkdir(parents=True, exist_ok=True)
multiqc_dir.mkdir(parents=True, exist_ok=True)

fastqs = []
with metadata.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        for key in ("fastq_1", "fastq_2"):
            value = row.get(key, "").strip()
            if value and not value.startswith("FILL_"):
                fastqs.append(value)

if not fastqs:
    raise SystemExit("No FASTQ paths found in metadata.")

subprocess.run(["fastqc", "--threads", threads, "--outdir", str(fastqc_dir), *fastqs], check=True)
subprocess.run(["multiqc", str(fastqc_dir), "--outdir", str(multiqc_dir)], check=True)
PY
