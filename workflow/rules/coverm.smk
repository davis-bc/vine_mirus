if READS_SAMPLES:
    rule coverm:
        input:
            viral_fa="results/viral_contigs.fa",
            r1=expand("results/samples/{sample}/qc/{sample}_R1.fastq.gz", sample=READS_SAMPLES),
            r2=expand("results/samples/{sample}/qc/{sample}_R2.fastq.gz", sample=READS_SAMPLES),
        output:
            abundance="results/intermediate/abundance/abundance_matrix.tsv",
        threads: 8
        resources:
            mem_mb=32000,
            runtime=120,
        params:
            metric=config["coverm_metric"],
            samples=READS_SAMPLES,
            test_mode=config.get("test_mode", False),
        conda:
            "../../envs/viral_analysis.yaml"
        log:
            "logs/coverm.log"
        script:
            "../scripts/run_coverm.py"
else:
    rule coverm:
        output:
            abundance="results/intermediate/abundance/abundance_matrix.tsv",
        log:
            "logs/coverm.log"
        shell:
            "mkdir -p results/intermediate/abundance && printf 'contig_id\\n' > {output.abundance} && echo '[no reads] skipped CoverM' > {log}"


