# Test Data

`data/` contains small synthetic paired-end FASTQ files for end-to-end pipeline testing.

## Generating test data

```bash
python test/generate_test_data.py
```

This creates:
- `data/sample1_R1.fastq.gz` / `data/sample1_R2.fastq.gz`
- `data/sample2_R1.fastq.gz` / `data/sample2_R2.fastq.gz`

Each sample has ~10,000 read pairs (150 bp), including ~1,000 reads derived from a
synthetic 50 kb phage-like genome to give the viral detection tools signal to work with.

## Running the test pipeline

```bash
# From repo root
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile \
          --config samples=config/samples.tsv
```

Expected runtime: < 10 minutes on a laptop (excluding conda environment setup).
