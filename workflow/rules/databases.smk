# databases.smk — automatic download/setup of all tool databases.
#
# Each rule produces a sentinel (.done) file or the actual DB artifact.
# Databases are stored under databases/ at the repo root and only downloaded once.
# Tool rules declare these sentinels as conditional inputs (skipped in test_mode).

# ── geNomad ───────────────────────────────────────────────────────────────────
rule download_genomad_db:
    """Download the geNomad database (~3 GB)."""
    output:
        done=touch("databases/genomad_db/.done"),
    params:
        outdir="databases/genomad_db",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=120,
    conda:
        "../../envs/genomad.yaml"
    log:
        "logs/databases/genomad_db.log"
    shell:
        "mkdir -p {params.outdir} && "
        "genomad download-database {params.outdir} 2> {log}"


# ── VirSorter2 ────────────────────────────────────────────────────────────────
rule download_virsorter2_db:
    """Download and set up the VirSorter2 database (~10 GB)."""
    output:
        done=touch("databases/virsorter2_db/.done"),
    params:
        outdir="databases/virsorter2_db",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=360,
    conda:
        "../../envs/virsorter2.yaml"
    log:
        "logs/databases/virsorter2_db.log"
    shell:
        "mkdir -p {params.outdir} && "
        "virsorter setup --db-dir {params.outdir} --jobs {threads} 2> {log}"


# ── CheckV ────────────────────────────────────────────────────────────────────
rule download_checkv_db:
    """Download the CheckV database (~3 GB)."""
    output:
        done=touch("databases/checkv_db/.done"),
    params:
        outdir="databases/checkv_db",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=120,
    conda:
        "../../envs/viral_analysis.yaml"
    log:
        "logs/databases/checkv_db.log"
    shell:
        "mkdir -p {params.outdir} && "
        "checkv download_database {params.outdir} 2> {log}"


# ── NCBI Taxdump ──────────────────────────────────────────────────────────────
rule download_taxdump:
    """Download NCBI taxonomy dump (nodes.dmp + names.dmp, ~70 MB compressed)."""
    output:
        done=touch("databases/taxdump/.done"),
    params:
        outdir="databases/taxdump",
        url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
    log:
        "logs/databases/taxdump.log"
    shell:
        """
        mkdir -p {params.outdir}
        curl -L {params.url} -o {params.outdir}/taxdump.tar.gz 2> {log}
        tar -xzf {params.outdir}/taxdump.tar.gz -C {params.outdir} 2>> {log}
        rm {params.outdir}/taxdump.tar.gz
        """


# ── DIAMOND viral RefSeq ──────────────────────────────────────────────────────
rule build_diamond_db:
    """Download NCBI viral RefSeq proteins and build a DIAMOND database (~2 GB)."""
    output:
        dmnd="databases/diamond/nr_viral.dmnd",
    params:
        outdir="databases/diamond",
        ftp_base="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=480,
    conda:
        "../../envs/viral_analysis.yaml"
    log:
        "logs/databases/diamond_db.log"
    shell:
        """
        mkdir -p {params.outdir}
        # Download all available viral RefSeq protein FASTA files
        for i in $(seq 1 100); do
            file="viral.${{i}}.protein.faa.gz"
            if curl -f -s -o "{params.outdir}/$file" \
                    "{params.ftp_base}/$file" 2>> {log}; then
                echo "Downloaded $file" >> {log}
            else
                rm -f "{params.outdir}/$file"
                break
            fi
        done
        # Build DIAMOND database
        cat {params.outdir}/viral.*.protein.faa.gz | \
            diamond makedb --in /dev/stdin \
            --db {params.outdir}/nr_viral \
            --threads {threads} 2>> {log}
        rm -f {params.outdir}/viral.*.protein.faa.gz
        """
