rule megahit:
    input:
        r1="results/samples/{sample}/qc/{sample}_R1.fastq.gz",
        r2="results/samples/{sample}/qc/{sample}_R2.fastq.gz",
    output:
        contigs="results/samples/{sample}/assembly/{sample}_contigs.fa",
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
    params:
        outdir="results/samples/{sample}/assembly/tmp",
        min_len=config["megahit_min_contig_len"],
    conda:
        "../../envs/qc_assembly.yaml"
    log:
        "logs/megahit/{sample}.log"
    shell:
        """
        rm -rf {params.outdir}
        megahit \
            -1 {input.r1} -2 {input.r2} \
            --min-contig-len {params.min_len} \
            -t {threads} \
            -m {resources.mem_mb}000000 \
            -o {params.outdir} \
            2> {log}
        cp {params.outdir}/final.contigs.fa {output.contigs}
        rm -rf {params.outdir}
        """
