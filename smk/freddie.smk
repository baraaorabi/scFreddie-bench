rule freddie:
    input:
        script="freddie.py",
        bam=f"{output_d}/{{sample}}/{{sample}}.sorted.bam",
        rname_to_celltypes=f"{output_d}/rname_to_celltypes/{{sample}}.tsv",
    output:
        isoforms=f"{output_d}/{{sample}}/freddie.isoforms.gtf",
    threads: 8
    resources:
        mem="16G",
        time=359,
    shell:
        "./{input.script}"
        " --rname-to-celltypes {input.rname_to_celltypes}"
        " --bam {input.bam}"
        " --output {output.isoforms}"
        " --threads {threads}"
