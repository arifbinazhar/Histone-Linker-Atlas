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
outdir = Path(config["paths"]["raw_fastq_dir"])
threads = str(config.get("tools", {}).get("threads", 8))
outdir.mkdir(parents=True, exist_ok=True)

with metadata.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        accession = row.get("sra_accession", "").strip()
        sample_id = row.get("sample_id", accession).strip()
        if not accession or accession.startswith("FILL_"):
            print(f"Skipping {sample_id}: fill sra_accession from GEO/SRA run table first.")
            continue
        print(f"Downloading {sample_id} ({accession})")
        subprocess.run(["prefetch", accession, "--output-directory", str(outdir / "sra")], check=True)
        subprocess.run(
            [
                "fasterq-dump",
                accession,
                "--split-files",
                "--threads",
                threads,
                "--outdir",
                str(outdir),
            ],
            check=True,
        )
PY
