rule evidence_score:
    input:
        genomad=expand("results/samples/{sample}/genomad/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv", sample=SAMPLES),
        virsorter2=expand("results/samples/{sample}/virsorter2/final-viral-score.tsv", sample=SAMPLES),
        pprmeta=expand("results/samples/{sample}/pprmeta/{sample}_pprmeta.csv", sample=SAMPLES),
    output:
        scores="results/intermediate/evidence/viral_scores.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        samples=SAMPLES,
        threshold=config["viral_score_threshold"],
        weights=config["tool_weights"],
    conda:
        "../../envs/scripts.yaml"
    log:
        "logs/evidence_score.log"
    script:
        "../scripts/score_evidence.py"
