rule fastp:
    input:
        r1=get_r1,
        r2=get_r2,
    output:
        r1="results/samples/{sample}/qc/{sample}_R1.fastq.gz",
        r2="results/samples/{sample}/qc/{sample}_R2.fastq.gz",
        json="results/samples/{sample}/qc/{sample}_fastp.json",
        html="results/samples/{sample}/qc/{sample}_fastp.html",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        min_len=config["fastp_min_length"],
        qual=config["fastp_qualified_quality"],
    conda:
        "../../envs/qc_assembly.yaml"
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --json {output.json} --html {output.html} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.qual} \
            --thread {threads} \
            2> {log}
        """
