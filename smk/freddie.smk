output_d = f"{config['outpath']}/freddie"


rule isoforms:
    input:
        bam=f"{output_d}/preprocess/{{sample}}.sorted.bam",
        rname_to_celltypes=f"{output_d}/preprocess/{{sample}}.rname_to_celltypes.tsv",
    output:
        isoforms=f"{output_d}/{{sample}}.isoforms.gtf",
    threads: 32
    params:
        script=config["freddie"]["exec"],
    resources:
        mem="32G",
        time=359,
    conda:
        "envs/freddie.yaml"
    shell:
        "{params.script}"
        " --rname-to-celltypes {input.rname_to_celltypes}"
        " --bam {input.bam}"
        " --output {output.isoforms}"
        " --threads {threads}"


rule rname_to_celltypes:
    input:
        lr_tsv=f"{output_d}/preprocess/{{sample}}.lr_matches.tsv.gz",
        truth=lambda wc: config["samples"][wc.sample]["truth"],
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.rname_to_celltypes.tsv",
    params:
        script=config["freddie"]["rname_to_celltypes"],
    conda:
        "envs/freddie.yaml"
    shell:
        "python {params.script}"
        " -lr {input.lr_tsv}"
        " -truth {input.truth}"
        " -o {output.tsv}"


rule scTagger_match:
    input:
        lr_tsv=f"{output_d}/preprocess/{{sample}}.lr_bc.tsv.gz",
        wl_tsv=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv.gz",
    output:
        lr_tsv=f"{output_d}/preprocess/{{sample}}.lr_matches.tsv.gz",
    threads: 32
    conda:
        "envs/sctagger.yaml"
    shell:
        "scTagger.py match_trie"
        " -lr {input.lr_tsv}"
        " -sr {input.wl_tsv}"
        " -o {output.lr_tsv}"
        " -t {threads}"


rule scTagger_extract_bc:
    input:
        tsv=f"{output_d}/preprocess/{{sample}}.lr_bc.tsv.gz",
        cb=lambda wc: config["samples"][wc.sample]["CB"],
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv.gz",
    conda:
        "envs/sctagger.yaml"
    shell:
        "scTagger.py extract_sr_bc_from_lr"
        " -i {input.tsv}"
        " -wl {input.cb}"
        " -o {output.tsv}"


rule scTagger_lr_seg:
    input:
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.lr_bc.tsv.gz",
    threads: 32
    conda:
        "envs/sctagger.yaml"
    shell:
        "scTagger.py extract_lr_bc"
        " -r {input.fastq}"
        " -o {output.tsv}"
        " -t {threads}"


rule minimap2:
    input:
        reads=lambda wc: config["samples"][wc.sample]["FASTQ"],
        genome=lambda wc: config["samples"][wc.sample]["DNA"],
    output:
        bam=f"{output_d}/preprocess/{{sample}}.sorted.bam",
        bai=f"{output_d}/preprocess/{{sample}}.sorted.bam.bai",
    threads: 32
    resources:
        mem="128G",
        time=24 * 60 - 1,
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -a -x splice -t {threads} {input.genome} {input.reads} | "
        "  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam - > {output.bam} && "
        "  samtools index {output.bam}"
