output_d = f"{config['outpath']}/scNanoGPS"


rule read_qc:
    input:
        script=config["scNanoGPS"]["read_qc"],
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
    output:
        o1=f"{output_d}/{{sample}}/first_tail.fastq.gz",
        o2=f"{output_d}/{{sample}}/last_tail.fastq.gz",
    params:
        d=f"{output_d}/{{sample}}",
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -i {input.fastq}"
        " -d {params.d}"


rule profile:
    input:
        script=config["scNanoGPS"]["profile"],
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
    output:
        png=f"{output_d}/{{sample}}/read_length.png",
    params:
        d=f"{output_d}/{{sample}}",
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -i {input.fastq}"
        " -d {params.d}"


rule scan:
    input:
        script=config["scNanoGPS"]["scan"],
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
    output:
        fastq=f"{output_d}/{{sample}}/processed.fastq.gz",
        cb=f"{output_d}/{{sample}}/barcode_list.tsv.gz",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -i {input.fastq}"
        " -d {params.d}"
        " -t {threads}"


rule assign:
    input:
        script=config["scNanoGPS"]["assign"],
        cb=f"{output_d}/{{sample}}/barcode_list.tsv.gz",
    output:
        cb_count=f"{output_d}/{{sample}}/CB_counting.tsv.gz",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -d {params.d}"
        " -t {threads}"


rule curate:
    input:
        script=config["scNanoGPS"]["curate"],
        fastq=f"{output_d}/{{sample}}/processed.fastq.gz",
        dna=lambda wc: config["samples"][wc.sample]["DNA"],
        cb_count=f"{output_d}/{{sample}}/CB_counting.tsv.gz",
        ref_genome=lambda wc: config["samples"][wc.sample]["DNA"],
    output:
        tmp_dir=directory(f"{output_d}/{{sample}}/tmp"),
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -d {params.d}"
        " --ref_genome={input.dna}"
        " --tmp_dir={output.tmp_dir}"
        " -t {threads}"


rule expression:
    input:
        script=config["scNanoGPS"]["expression"],
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
        tmp_dir=f"{output_d}/{{sample}}/tmp",
    output:
        filtered_barcode_list=f"{output_d}/{{sample}}/filtered_barcode_list.txt",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -d {params.d}"
        " --gtf={input.gtf}"
        " --tmp_dir={input.tmp_dir}"
        " -t {threads}"


rule liqa_refgene:
    input:
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
    output:
        liqa_refgene=f"{output_d}/{{sample}}/liqa.refgene",
    shell:
        "liqa -task refgene"
        " -format gtf"
        " -ref {input.gtf}"
        " -out {output.liqa_refgene}"


rule isoforms:
    input:
        script=config["scNanoGPS"]["isoforms"],
        filtered_barcode_list=f"{output_d}/{{sample}}/filtered_barcode_list.txt",
        liqa_refgene=f"{output_d}/{{sample}}/liqa.refgene",
        tmp_dir=f"{output_d}/{{sample}}/tmp",
    output:
        matrix_isoform=f"{output_d}/{{sample}}/matrix_isoform.tsv",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " --liqa_ref {input.liqa_refgene}"
        " --tmp_dir={input.tmp_dir}"
        " -d {params.d}"
        " -t {threads}"
