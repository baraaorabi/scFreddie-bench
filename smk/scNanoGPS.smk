output_d = f"{config['outpath']}/scNanoGPS"


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
        filtered_barcode_list=f"{output_d}/{{sample}}/filtered_barcode_list.txt",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -d {params.d}"
        " --ref_genome={input.dna}"
        " --tmp_dir={params.d}/tmp"
        " --skip_curation 1"
        " -t {threads}"


rule expression:
    input:
        script=config["scNanoGPS"]["expression"],
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
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
        " -t {threads}"


rule report:
    input:
        script=config["scNanoGPS"]["report"],
        filtered_barcode_list=f"{output_d}/{{sample}}/filtered_barcode_list.txt",
        liqa_ref=lambda wc: config["samples"][wc.sample]["DNA"],
    output:
        matrix_isoform=f"{output_d}/{{sample}}/matrix_isoform.tsv",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        "python {input.script}"
        " -d {params.d}"
        " -t {threads}"
