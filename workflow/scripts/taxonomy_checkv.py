"""
taxonomy_checkv.py — extract viral contigs, run DIAMOND + taxopy for taxonomy,
then run CheckV for completeness. Called as a Snakemake script.

Combines taxonomy assignment and CheckV into one rule since CheckV is lightweight.
"""

import sys
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import taxopy

# ── Setup ─────────────────────────────────────────────────────────────────────
log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

scores_path  = snakemake.input.scores
contig_paths = snakemake.input.contigs

diamond_db   = snakemake.params.diamond_db
taxdump_dir  = snakemake.params.taxdump
checkv_db    = snakemake.params.checkv_db
evalue       = snakemake.params.evalue
max_target   = snakemake.params.max_target
outfmt       = snakemake.params.outfmt
checkv_outdir = snakemake.params.checkv_outdir
min_len      = snakemake.params.min_len
threads      = snakemake.threads

viral_fa_out    = snakemake.output.viral_fa
taxonomy_out    = snakemake.output.taxonomy
checkv_out      = snakemake.output.checkv

# ── Step 1: Extract passing viral contigs ─────────────────────────────────────
scores = pd.read_csv(scores_path, sep="\t")
passing = set(scores.loc[scores["PASS"], "contig_id"].tolist())

print(f"Extracting {len(passing)} viral contigs from {len(contig_paths)} assemblies...")

Path(viral_fa_out).parent.mkdir(parents=True, exist_ok=True)

contig_to_sample = {}
with open(viral_fa_out, "w") as out_fa:
    for contig_path in contig_paths:
        sample = Path(contig_path).parent.parent.name  # results/{sample}/megahit/
        for rec in SeqIO.parse(contig_path, "fasta"):
            if rec.id in passing and len(rec.seq) >= min_len:
                SeqIO.write(rec, out_fa, "fasta")
                contig_to_sample[rec.id] = sample

print(f"Wrote {len(contig_to_sample)} sequences to {viral_fa_out}")

# ── Step 2: DIAMOND blastp for taxonomy ───────────────────────────────────────
Path(taxonomy_out).parent.mkdir(parents=True, exist_ok=True)
diamond_out = Path(taxonomy_out).parent / "diamond_hits.tsv"

print("Running DIAMOND blastp...")
subprocess.run(
    [
        "diamond", "blastp",
        "--db", diamond_db,
        "--query", viral_fa_out,
        "--out", str(diamond_out),
        "--outfmt", *outfmt.split(),
        "--evalue", str(evalue),
        "--max-target-seqs", str(max_target),
        "--threads", str(threads),
        "--more-sensitive",
    ],
    check=True, stderr=log,
)

# ── Step 3: Resolve lineages with taxopy ──────────────────────────────────────
print("Loading NCBI taxdump...")
taxdb = taxopy.TaxDb(
    nodes_dmp=str(Path(taxdump_dir) / "nodes.dmp"),
    names_dmp=str(Path(taxdump_dir) / "names.dmp"),
)

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

def taxid_to_lineage(taxid):
    try:
        taxon = taxopy.Taxon(int(taxid), taxdb)
        lineage = {r: "" for r in RANKS}
        for name, rank in zip(taxon.name_lineage, taxon.rank_lineage):
            if rank in lineage:
                lineage[rank] = name
        return ";".join(lineage[r] for r in RANKS)
    except Exception:
        return ";;;;;;;"

hits = pd.read_csv(
    diamond_out, sep="\t", header=None,
    names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "sscinames"],
)
# Best hit per query (first row after sort by bitscore)
best = hits.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
best["taxonomy"] = best["staxids"].apply(lambda x: taxid_to_lineage(str(x).split(";")[0]))
best["sample"] = best["qseqid"].map(contig_to_sample)

tax_df = best[["qseqid", "sample", "taxonomy", "pident", "evalue", "bitscore"]].rename(
    columns={"qseqid": "contig_id"}
)
# Add contigs with no DIAMOND hit (score too low)
all_contigs = pd.DataFrame(
    {"contig_id": list(contig_to_sample.keys()), "sample": list(contig_to_sample.values())}
)
tax_df = all_contigs.merge(tax_df, on=["contig_id", "sample"], how="left")
tax_df["taxonomy"] = tax_df["taxonomy"].fillna(";;;;;;;")
tax_df.to_csv(taxonomy_out, sep="\t", index=False)
print(f"Taxonomy written to {taxonomy_out}")

# ── Step 4: CheckV ────────────────────────────────────────────────────────────
print("Running CheckV...")
Path(checkv_outdir).mkdir(parents=True, exist_ok=True)
subprocess.run(
    [
        "checkv", "end_to_end",
        viral_fa_out, checkv_outdir,
        "-d", checkv_db,
        "-t", str(threads),
    ],
    check=True, stderr=log,
)
print(f"CheckV complete. Summary at {checkv_out}")
log.close()
