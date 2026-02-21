PPRMETA_DIR = "databases/pprmeta"
PPRMETA_SCRIPT = f"{PPRMETA_DIR}/PPR_Meta.py"
TEST_MODE = config.get("test_mode", False)
_SCRIPTS = workflow.basedir + "/scripts"

_PPRMETA_CMD = (
    f"python {_SCRIPTS}/mock_pprmeta.py "
    "{input.contigs} {output.results} {params.min_len} "
    "2> {log}"
    if TEST_MODE else
    "mkdir -p {params.outdir} && "
    f"python {_SCRIPTS}/filter_fasta_by_length.py "
    "{input.contigs} {params.outdir}/filtered.fa {params.min_len} && "
    "python3 {input.script} {params.outdir}/filtered.fa {output.results} "
    "2> {log}"
)


rule setup_pprmeta:
    """Download PPR-Meta v1.1 from GitHub into databases/pprmeta/ if not present."""
    output:
        script=PPRMETA_SCRIPT,
    params:
        url="https://github.com/zhenchengfang/PPR-Meta/releases/download/v1.1/PPR_Meta_v_1_1.zip",
        outdir=PPRMETA_DIR,
    log:
        "logs/setup_pprmeta.log"
    shell:
        """
        mkdir -p {params.outdir}
        curl -L {params.url} -o {params.outdir}/PPR_Meta.zip 2> {log}
        unzip -o {params.outdir}/PPR_Meta.zip -d {params.outdir} >> {log} 2>&1
        find {params.outdir} -name "PPR_Meta.py" ! -path "{output.script}" \
            -exec mv {{}} {output.script} \\; 2>> {log} || true
        """


def pprmeta_inputs(wildcards):
    """In test_mode skip the setup_pprmeta download; script input not needed."""
    base = {"contigs": get_contigs_path(wildcards.sample)}
    if not TEST_MODE:
        base["script"] = PPRMETA_SCRIPT
    return base


rule pprmeta:
    input:
        unpack(pprmeta_inputs),
    output:
        results="results/samples/{sample}/pprmeta/{sample}_pprmeta.csv",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        outdir="results/samples/{sample}/pprmeta",
        min_len=config["min_contig_length"],
    conda:
        "../../envs/pprmeta.yaml"
    log:
        "logs/pprmeta/{sample}.log"
    shell:
        _PPRMETA_CMD
