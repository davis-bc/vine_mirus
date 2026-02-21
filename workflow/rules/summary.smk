rule summary:
    input:
        scores="results/intermediate/evidence/viral_scores.tsv",
        taxonomy="results/intermediate/taxonomy/viral_taxonomy.tsv",
        checkv="results/intermediate/checkv/quality_summary.tsv",
        abundance="results/intermediate/abundance/abundance_matrix.tsv",
    output:
        tsv="results/viral_contigs_abundance.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        samples=READS_SAMPLES,
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/summary.log"
    script:
        "../scripts/summarize.py"
