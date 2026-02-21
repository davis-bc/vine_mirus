"""
mock_genomad.py — generate synthetic geNomad virus_summary.tsv output for test_mode.
Called via shell: python {workflow}/scripts/mock_genomad.py <contigs> <output> <min_len>
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
        "seq_name": rec.id,
        "length": len(rec.seq),
        "topology": "linear",
        "coordinates": f"1-{len(rec.seq)}",
        "n_genes": random.randint(1, 10),
        "genetic_code": 11,
        "virus_score": score,
        "fdr": round(random.uniform(0, 0.5), 4),
        "n_hallmarks": random.randint(0, 3),
        "marker_enrichment": round(random.uniform(0, 5), 4),
        "taxonomy": "Viruses;Duplodnaviria;;Caudoviricetes;;;" if score > 0.7 else "",
    })

pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
print(f"[mock] geNomad: wrote {len(rows)} rows to {args.output}")
