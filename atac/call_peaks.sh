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
filtered_dir = Path(config["paths"].get("filtered_dir", "results/atac/GSE236992/filtered"))
peak_dir = Path(config["paths"]["peak_dir"])
qvalue = str(config.get("tools", {}).get("macs2_qvalue", 0.05))
peak_dir.mkdir(parents=True, exist_ok=True)

peak_files = []
with metadata.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        sample = row["sample_id"]
        bam = filtered_dir / f"{sample}.filtered.bam"
        out_prefix = peak_dir / sample
        subprocess.run(
            [
                "macs2",
                "callpeak",
                "-t",
                str(bam),
                "-f",
                "BAMPE",
                "-g",
                "hs",
                "-n",
                sample,
                "--outdir",
                str(peak_dir),
                "--nomodel",
                "--shift",
                "-100",
                "--extsize",
                "200",
                "-q",
                qvalue,
            ],
            check=True,
        )
        peak_files.append(str(out_prefix) + "_peaks.narrowPeak")

consensus = peak_dir / "consensus_peaks.bed"
cmd = f"cat {' '.join(peak_files)} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > {consensus}"
subprocess.run(cmd, shell=True, check=True)
PY
