"""
mock_virsorter2.py — generate synthetic VirSorter2 final-viral-score.tsv for test_mode.
Called via shell: python {workflow}/scripts/mock_virsorter2.py <contigs> <output> <min_len>
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
    score = round(random.uniform(0.0, 1.0), 4)
    rows.append({
        "seqname": f"{rec.id}||full",
        "dsDNAphage": round(score * 0.8, 4),
        "ssDNA": round(score * 0.1, 4),
        "NCLDV": round(score * 0.05, 4),
        "RNA": round(score * 0.03, 4),
        "lavidaviridae": round(score * 0.02, 4),
        "max_score": score,
        "max_score_group": "dsDNAphage",
        "length": len(rec.seq),
        "hallmark": random.randint(0, 2),
        "viral": "Yes" if score >= 0.5 else "No",
        "cellular": "No" if score >= 0.5 else "Yes",
    })

pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
print(f"[mock] VirSorter2: wrote {len(rows)} rows to {args.output}")
