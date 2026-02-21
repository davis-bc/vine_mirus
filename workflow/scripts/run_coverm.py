"""
run_coverm.py — compute per-contig read abundance using CoverM (or mock in test_mode).
Called via Snakemake script: directive with viral_analysis conda env.
"""

import random
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

random.seed(42)

log_fh = open(snakemake.log[0], "w")
sys.stdout = log_fh
sys.stderr = log_fh

params    = snakemake.params
test_mode = params.test_mode
threads   = snakemake.threads

out_path = snakemake.output.abundance
Path(out_path).parent.mkdir(parents=True, exist_ok=True)

if test_mode:
    print("[mock] CoverM")
    contigs = list(SeqIO.parse(snakemake.input.viral_fa, "fasta"))
    rows = []
    for rec in contigs:
        row = {"Contig": rec.id}
        for s in params.samples:
            row[f"{s} RPKM"] = round(random.uniform(0.0, 1000.0), 2)
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    print(f"[mock] CoverM: {len(rows)} contigs x {len(params.samples)} samples")
else:
    from snakemake.shell import shell
    r1 = snakemake.input.r1
    r2 = snakemake.input.r2
    coupled = " ".join(f"{a} {b}" for a, b in zip(r1, r2))
    shell(
        "coverm contig "
        "--reference {snakemake.input.viral_fa} "
        f"--coupled {coupled} "
        "--methods {params.metric} "
        "--threads {threads} "
        "--output-file {out_path} "
        "2>> {snakemake.log[0]}"
    )

log_fh.close()
