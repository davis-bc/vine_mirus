_SCRIPTS = workflow.basedir + "/scripts"

_GENOMAD_CMD = (
    f"python {_SCRIPTS}/mock_genomad.py "
    "{input.contigs} {output.summary} {params.min_len} "
    "2> {log}"
    if config.get("test_mode") else
    "genomad end-to-end "
    "--min-score 0.0 "
    "--min-length {params.min_len} "
    "--threads {threads} "
    "--cleanup "
    "{input.contigs} {params.outdir} {params.db} "
    "2> {log}"
)

rule genomad:
    input:
        contigs=get_contigs,
        db=[] if config.get("test_mode") else ["databases/genomad_db/.done"],
    output:
        summary="results/samples/{sample}/genomad/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        outdir="results/samples/{sample}/genomad",
        db=config["genomad_db"],
        min_len=config["min_contig_length"],
    conda:
        "../../envs/genomad.yaml"
    log:
        "logs/genomad/{sample}.log"
    shell:
        _GENOMAD_CMD
