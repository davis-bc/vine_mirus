"""
score_evidence.py — fuse viral predictions from geNomad, VirSorter2, and PPR-Meta
into a single probability score per contig. Called as a Snakemake script.

Inputs  (snakemake.input):
  genomad   : list of geNomad virus_summary.tsv paths (one per sample)
  virsorter2: list of VirSorter2 final-viral-score.tsv paths (one per sample)
  pprmeta   : list of PPR-Meta CSV output paths (one per sample)

Output (snakemake.output.scores):
  TSV with columns: contig_id, sample, genomad_score, virsorter2_score,
                    pprmeta_score, viral_score, PASS
"""

import sys
import pandas as pd
from pathlib import Path

samples      = snakemake.params.samples
threshold    = snakemake.params.threshold
weights      = snakemake.params.weights
log_path     = snakemake.log[0]

sys.stderr = open(log_path, "w")

w_gn = weights["genomad"]
w_vs = weights["virsorter2"]
w_pp = weights["pprmeta"]


def load_genomad(path, sample):
    try:
        df = pd.read_csv(path, sep="\t")
        # geNomad column is 'virus_score' or 'score'
        score_col = "virus_score" if "virus_score" in df.columns else "score"
        df = df[["seq_name", score_col]].rename(columns={"seq_name": "contig_id", score_col: "genomad_score"})
        df["sample"] = sample
        return df
    except Exception as e:
        print(f"WARNING: could not parse geNomad output for {sample}: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["contig_id", "sample", "genomad_score"])


def load_virsorter2(path, sample):
    try:
        df = pd.read_csv(path, sep="\t")
        # strip VS2 suffix (||full, ||partial, etc.)
        df["contig_id"] = df["seqname"].str.replace(r"\|\|.*$", "", regex=True)
        df = df[["contig_id", "max_score"]].rename(columns={"max_score": "virsorter2_score"})
        df["sample"] = sample
        return df
    except Exception as e:
        print(f"WARNING: could not parse VirSorter2 output for {sample}: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["contig_id", "sample", "virsorter2_score"])


def load_pprmeta(path, sample):
    try:
        df = pd.read_csv(path)
        # PPR-Meta columns: Header, Phage_score, Chromosome_score, Plasmid_score, Possible_source
        df["contig_id"] = df["Header"].str.split(r"\s+").str[0]
        df = df[["contig_id", "Phage_score"]].rename(columns={"Phage_score": "pprmeta_score"})
        df["sample"] = sample
        return df
    except Exception as e:
        print(f"WARNING: could not parse PPR-Meta output for {sample}: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["contig_id", "sample", "pprmeta_score"])


genomad_dfs    = [load_genomad(p, s) for p, s in zip(snakemake.input.genomad,    samples)]
virsorter2_dfs = [load_virsorter2(p, s) for p, s in zip(snakemake.input.virsorter2, samples)]
pprmeta_dfs    = [load_pprmeta(p, s) for p, s in zip(snakemake.input.pprmeta,    samples)]

gn = pd.concat(genomad_dfs,    ignore_index=True)
vs = pd.concat(virsorter2_dfs, ignore_index=True)
pp = pd.concat(pprmeta_dfs,    ignore_index=True)

# Merge on contig_id + sample; outer join so contigs detected by any tool are kept
merged = gn.merge(vs, on=["contig_id", "sample"], how="outer") \
           .merge(pp, on=["contig_id", "sample"], how="outer")

merged["genomad_score"]    = merged["genomad_score"].fillna(0.0)
merged["virsorter2_score"] = merged["virsorter2_score"].fillna(0.0)
merged["pprmeta_score"]    = merged["pprmeta_score"].fillna(0.0)

merged["viral_score"] = (
    w_gn * merged["genomad_score"] +
    w_vs * merged["virsorter2_score"] +
    w_pp * merged["pprmeta_score"]
)

merged["PASS"] = merged["viral_score"] >= threshold

# Reorder columns: sample | contig_id | per-tool scores | viral_score | PASS
merged = merged[["sample", "contig_id", "genomad_score", "virsorter2_score", "pprmeta_score", "viral_score", "PASS"]]

out = snakemake.output.scores
Path(out).parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(out, sep="\t", index=False)

n_pass = merged["PASS"].sum()
print(f"Total contigs scored: {len(merged)}; passing threshold ({threshold}): {n_pass}", file=sys.stderr)
