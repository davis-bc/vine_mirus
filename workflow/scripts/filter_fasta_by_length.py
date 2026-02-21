"""filter_fasta_by_length.py — filter a FASTA to records >= min_len bp."""

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("min_len", type=int)
args = parser.parse_args()

records = [r for r in SeqIO.parse(args.input, "fasta") if len(r.seq) >= args.min_len]
SeqIO.write(records, args.output, "fasta")
print(f"Filtered {len(records)} records >= {args.min_len} bp to {args.output}")
