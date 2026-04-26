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
trimmed_dir = Path(config["paths"]["trimmed_dir"])
aligned_dir = Path(config["paths"]["aligned_dir"])
index = config["reference"]["bowtie2_index"]
threads = str(config.get("tools", {}).get("threads", 8))
aligned_dir.mkdir(parents=True, exist_ok=True)

with metadata.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        sample = row["sample_id"]
        fq1 = trimmed_dir / f"{sample}.trimmed.R1.fastq.gz"
        fq2 = trimmed_dir / f"{sample}.trimmed.R2.fastq.gz"
        bam = aligned_dir / f"{sample}.sorted.bam"
        log = aligned_dir / f"{sample}.bowtie2.log"
        cmd = (
            f"bowtie2 -x {index} -1 {fq1} -2 {fq2} --very-sensitive -X 2000 -p {threads} "
            f"2> {log} | samtools sort -@ {threads} -o {bam} -"
        )
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
PY
