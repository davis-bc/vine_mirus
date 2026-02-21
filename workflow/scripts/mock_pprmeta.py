"""
mock_pprmeta.py — generate synthetic PPR-Meta CSV output for test_mode.
Called via shell: python {workflow}/scripts/mock_pprmeta.py <contigs> <output> <min_len>
"""

import argparse
import random
from pathlib import Path

import pandas as pd
from Bio import SeqIO

random.seed(42)

parser = argparse.ArgumentParser()
parser.add_argument("contigs")
parser.add_argument("output")
parser.add_argument("min_len", type=int)
args = parser.parse_args()

contigs = [r for r in SeqIO.parse(args.contigs, "fasta") if len(r.seq) >= args.min_len]
Path(args.output).parent.mkdir(parents=True, exist_ok=True)

rows = []
for rec in contigs:
    phage = round(random.uniform(0.0, 1.0), 4)
    chrom = round(random.uniform(0.0, 1 - phage), 4)
    plas  = round(1 - phage - chrom, 4)
    rows.append({
        "Header": rec.id,
        "Phage_score": phage,
        "Chromosome_score": chrom,
        "Plasmid_score": plas,
        "Possible_source": "phage" if phage >= 0.5 else ("chromosome" if chrom >= 0.5 else "plasmid"),
    })

pd.DataFrame(rows).to_csv(args.output, index=False)
print(f"[mock] PPR-Meta: wrote {len(rows)} rows to {args.output}")
