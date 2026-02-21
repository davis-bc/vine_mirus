_SCRIPTS = workflow.basedir + "/scripts"

_VIRSORTER2_CMD = (
    f"python {_SCRIPTS}/mock_virsorter2.py "
    "{input.contigs} {output.scores} {params.min_len} "
    "2> {log}"
    if config.get("test_mode") else
    "virsorter run "
    "-i {input.contigs} "
    "-w {params.outdir} "
    "--include-groups {params.groups} "
    "--min-score {params.min_score} "
    "--min-length {params.min_len} "
    "-j {threads} "
    "--db-dir {params.db} "
    "all "
    "2> {log}"
)

rule virsorter2:
    input:
        contigs=get_contigs,
        db=[] if config.get("test_mode") else ["databases/virsorter2_db/.done"],
    output:
        scores="results/samples/{sample}/virsorter2/final-viral-score.tsv",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=180,
    params:
        outdir="results/samples/{sample}/virsorter2",
        db=config["virsorter2_db"],
        groups=",".join(config["virsorter2_groups"]),
        min_score=config["virsorter2_min_score"],
        min_len=config["min_contig_length"],
    conda:
        "../../envs/virsorter2.yaml"
    log:
        "logs/virsorter2/{sample}.log"
    shell:
        _VIRSORTER2_CMD
