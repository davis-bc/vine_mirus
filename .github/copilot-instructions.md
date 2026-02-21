# Copilot Instructions — viral_discovery

## Project Overview

A Snakemake bioinformatics pipeline for discovering viral sequences in metagenomic samples. Takes paired-end FASTQ files as input and produces a viral contig abundance matrix with taxonomy.

## Pipeline Architecture

The workflow is a linear DAG of stages. Each stage's outputs feed the next:

```
raw FASTQs (per sample)
  └─► fastp          → QC-filtered reads (per sample)
        └─► MEGAHIT  → assembled contigs (per sample)
              └─► [geNomad, VirSorter2, PPR-Meta]  → per-tool viral predictions (parallel)
                    └─► evidence scoring            → viral contig list with probability scores
                          └─► taxonomy assignment   → taxonomic lineage per contig
                                └─► CheckV          → completeness / quality assessment
                                      └─► CoverM    → RPKM abundance matrix (reads → contigs)
                                            └─► final summary TSV
```

**Multi-tool evidence fusion**: geNomad, VirSorter2, and PPR-Meta run in parallel on the same contigs. A scoring step combines their outputs into a single probability score per contig. Only contigs flagged by at least one tool proceed to downstream steps (configurable threshold in `config/config.yaml`).

## Repository Layout

```
workflow/
  Snakefile            # top-level rule (rule all targets)
  rules/               # one .smk file per pipeline stage
    fastp.smk
    megahit.smk
    genomad.smk
    virsorter2.smk
    pprmeta.smk
    evidence_score.smk
    taxonomy.smk
    checkv.smk
    coverm.smk
    summary.smk
  scripts/             # Python helper scripts called by rules
    score_evidence.py  # fuses multi-tool predictions into probability scores
    summarize.py       # builds final TSV
config/
  config.yaml          # all tuneable parameters (thresholds, threads, paths)
  samples.tsv          # sample name → R1/R2 FASTQ paths
resources/             # reference databases (not committed; paths in config.yaml)
results/               # all output; gitignored
  {sample}/
    fastp/
    megahit/
    genomad/
    virsorter2/
    pprmeta/
    evidence/
    taxonomy/
    checkv/
    coverm/
logs/                  # per-rule log files; gitignored
test/
  data/                # small synthetic FASTQs for end-to-end testing
  README.md            # how to generate/regenerate test data
envs/                  # Conda environment YAML files (one per tool where needed)
profiles/
  slurm/
    config.yaml        # snakemake-executor-plugin-slurm profile
```

## Snakemake Conventions

- **Run on cluster (SLURM)**: `snakemake --profile profiles/slurm`
- **Run locally (dev/test)**: `snakemake --cores all --use-conda`
- **Dry run**: `snakemake --profile profiles/slurm -n`
- **Single rule (local testing)**: `snakemake results/{sample}/fastp/{sample}_R1.fastq.gz --cores 4 --use-conda`
- **Run on test data locally**: `snakemake --cores 4 --use-conda --config samples=test/samples.tsv`
- **Force re-run a stage**: `snakemake --profile profiles/slurm --forcerun megahit`
- **DAG visualization**: `snakemake --dag | dot -Tpdf > dag.pdf`

All rules specify `threads:`, `resources:` (with `mem_mb` and `runtime`), and `log:` directives. Shell commands redirect stderr to `{log}`.

## Configuration

All thresholds and resource parameters live in `config/config.yaml`. Key entries:
- `min_contig_length` — minimum contig length passed to viral detection tools (default: 1000 bp)
- `viral_score_threshold` — minimum combined probability to retain a contig (default: 0.5)
- `coverm_metric` — abundance metric used (default: `rpkm`)
- Database paths for geNomad, VirSorter2, CheckV, and NCBI taxonomy (`taxdump/`)

`samples.tsv` columns: `sample`, `r1`, `r2` (tab-separated, one row per sample).

## Evidence Scoring Logic

`workflow/scripts/score_evidence.py` merges predictions from all three tools. Each tool contributes a score:
- **geNomad**: uses `score` column from `{sample}_summary/` output
- **VirSorter2**: uses `max_score` from `final-viral-score.tsv`
- **PPR-Meta**: uses `Possible_source` probability for phage/virus categories

Scores are averaged (equal weight by default; weights configurable in `config.yaml`). A contig passes if its combined score ≥ `viral_score_threshold`.

## Taxonomy

Use **NCBI taxonomy** throughout. The `taxonomy.smk` rule runs `diamond blastp` against the NCBI nr database (or a viral subset), then resolves lineages using the NCBI taxdump (`nodes.dmp`, `names.dmp`) via the `taxopy` Python library. geNomad also outputs NCBI taxids directly — prefer those when available and fall back to BLAST-based assignment. Lineage strings are formatted as semicolon-delimited NCBI ranks: `superkingdom;phylum;class;order;family;genus;species`.

## SLURM Execution

Use the `snakemake-executor-plugin-slurm` executor. The SLURM profile lives in `profiles/slurm/`:

```
profiles/slurm/
  config.yaml          # executor + default resources
```

**`profiles/slurm/config.yaml`:**
```yaml
executor: slurm
default-resources:
  slurm_partition: "main"
  mem_mb: 8000
  runtime: 60          # minutes
  tasks: 1
jobs: 200
latency-wait: 60
use-conda: true
```

**Run on the cluster:**
```bash
snakemake --profile profiles/slurm
```

**Dry run against SLURM:**
```bash
snakemake --profile profiles/slurm -n
```

### Per-rule SLURM resources

Override defaults in each rule's `resources:` block using `snakemake-executor-plugin-slurm` keys. Every rule must declare `mem_mb`, `runtime`, and `threads`. GPU rules additionally set `slurm_extra`.

| Rule | threads | mem_mb | runtime | notes |
|---|---|---|---|---|
| `fastp` | 8 | 8 000 | 60 | per sample |
| `megahit` | 16 | 64 000 | 240 | per sample; memory-intensive |
| `genomad` | 8 | 32 000 | 120 | per sample |
| `virsorter2` | 8 | 32 000 | 180 | per sample |
| `pprmeta` | 4 | 16 000 | 60 | per sample |
| `evidence_score` | 1 | 4 000 | 30 | aggregate step |
| `taxonomy` | 8 | 32 000 | 120 | DIAMOND + taxopy |
| `checkv` | 4 | 16 000 | 60 | per sample |
| `coverm` | 8 | 32 000 | 120 | per sample; includes BWA mem |
| `summary` | 1 | 4 000 | 30 | aggregate step |

Example rule skeleton:
```python
rule megahit:
    input: ...
    output: ...
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
        slurm_partition="main",
    conda: "../envs/megahit.yaml"
    log: "logs/megahit/{sample}.log"
    shell:
        "megahit -1 {input.r1} -2 {input.r2} "
        "--min-contig-len {params.min_len} "
        "-t {threads} --memory {resources.mem_mb}000000 "
        "-o {params.outdir} 2> {log}"
```

For jobs needing a high-memory or GPU partition, set `slurm_partition` explicitly in the rule's `resources:` block (e.g., `slurm_partition="highmem"`). Do not hard-code partition names in the profile default so they can be overridden per rule.

## Workflow Design Principles

**Minimize rules and Conda environments** to reduce setup overhead and Conda solve time.

### Rule consolidation
Combine steps that always run sequentially on the same input into a single rule rather than chaining many small rules. Prefer fewer, more self-contained rules over fine-grained ones. Only split into separate rules when a step is reused by multiple downstream rules, or when dramatically different resource requirements make a combined rule wasteful.

Preferred groupings:
- `fastp` + `megahit` can each remain single rules (different resource profiles), but intermediate FASTQ cleanup can happen inside the megahit rule rather than as a separate step.
- The three viral prediction tools (`genomad`, `virsorter2`, `pprmeta`) stay separate because they run in parallel and have incompatible dependencies, but their output-parsing logic should be absorbed into `evidence_score` rather than added as extra rules.
- `taxonomy` + `checkv` may run sequentially in one rule if resource profiles are compatible, since CheckV is lightweight.
- `coverm` encapsulates read mapping (BWA/minimap2) + coverage calculation internally — do not add a separate alignment rule.

### Conda environment consolidation
Target **4–5 environments total**. Group tools by dependency compatibility:

| Environment file | Tools |
|---|---|
| `envs/qc_assembly.yaml` | fastp, MEGAHIT |
| `envs/genomad.yaml` | geNomad (has heavy isolated deps) |
| `envs/virsorter2.yaml` | VirSorter2 (scikit-learn, imbalanced-learn) |
| `envs/viral_analysis.yaml` | PPR-Meta, CheckV, taxopy, DIAMOND, CoverM |
| `envs/scripts.yaml` | pandas, biopython — for `score_evidence.py` and `summarize.py` |

Before creating a new environment, check whether the tool can be added to an existing one. Only create a new `envs/*.yaml` when a genuine dependency conflict exists.

## Conda Environments

Environments are defined in `envs/` and referenced via `conda:` directives in rules. Do not install tools into the base environment.

## Test Dataset

`test/data/` contains small synthetic paired-end FASTQs (~10k reads per sample, 2 samples). Generate/regenerate with:
```bash
python test/generate_test_data.py
```
The test run should complete end-to-end in under 10 minutes on a laptop. Use it to validate any pipeline changes before touching real data.

## Output Format

The final output is `results/summary/viral_contigs_abundance.tsv` with columns:
`contig_id`, `sample_of_origin`, `length`, `viral_score`, `taxonomy` (semicolon-delimited lineage), `checkv_quality`, then one RPKM column per sample.
