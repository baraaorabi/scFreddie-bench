import format_time

output_d = f"{config['outpath']}/scNanoGPS"


rule make_time:
    output:
        time=f"{output_d}/time.tsv",
    run:
        with open(output.time, "w+") as fout:
            print(format_time.header_str, file=fout)


rule read_qc:
    input:
        script=config["scNanoGPS"]["read_qc"],
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        o1=f"{output_d}/{{sample}}/first_tail.fastq.gz",
        o2=f"{output_d}/{{sample}}/last_tail.fastq.gz",
    params:
        d=f"{output_d}/{{sample}}",
    conda:
        "envs/scnanogps.yaml"
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_read_qc')}"
        "python {input.script}"
        " -i {input.fastq}"
        " -d {params.d}"


rule profile:
    input:
        script=config["scNanoGPS"]["profile"],
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        png=f"{output_d}/{{sample}}/read_length.png",
    params:
        d=f"{output_d}/{{sample}}",
    conda:
        "envs/scnanogps.yaml"
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_profile')}"
        "python {input.script}"
        " -i {input.fastq}"
        " -d {params.d}"


rule scan:
    input:
        script=config["scNanoGPS"]["scan"],
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        fastq=f"{output_d}/{{sample}}/processed.fastq.gz",
        cb=f"{output_d}/{{sample}}/barcode_list.tsv.gz",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_scan')}"
        "python {input.script}"
        " -i {input.fastq}"
        " -d {params.d}"
        " -t {threads}"


rule assign:
    input:
        script=config["scNanoGPS"]["assign"],
        cb=f"{output_d}/{{sample}}/barcode_list.tsv.gz",
        time=ancient(f"{output_d}/time.tsv"),
    output:
        cb_count=f"{output_d}/{{sample}}/CB_counting.tsv.gz",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_assign')}"
        "python {input.script}"
        " -d {params.d}"
        " -t {threads}"


rule minimap2_index:
    input:
        fasta=lambda wc: config["samples"][wc.sample]["DNA"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        index=f"{output_d}/{{sample}}/ref_genome.mmi",
    threads: 32
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_minimap2_index')}"
        "minimap2 -t 32 -x map-ont -d {output.index} {input.fasta}"


rule curate:
    input:
        script=config["scNanoGPS"]["curate"],
        fastq=f"{output_d}/{{sample}}/processed.fastq.gz",
        dna=lambda wc: config["samples"][wc.sample]["DNA"],
        cb_count=f"{output_d}/{{sample}}/CB_counting.tsv.gz",
        ref_genome=lambda wc: config["samples"][wc.sample]["DNA"],
        index=f"{output_d}/{{sample}}/ref_genome.mmi",
        time=ancient(f"{output_d}/time.tsv"),
    output:
        tmp_dir=directory(f"{output_d}/{{sample}}/tmp"),
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    resources:
        mem="128G",
        time=72 * 60 - 1,
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_curate')}"
        "python {input.script}"
        " -d {params.d}"
        " --ref_genome {input.dna}"
        " --idx_genome {input.index}"
        " --tmp_dir {output.tmp_dir}"
        " -t {threads}"


rule expression:
    input:
        script=config["scNanoGPS"]["expression"],
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
        tmp_dir=f"{output_d}/{{sample}}/tmp",
        time=ancient(f"{output_d}/time.tsv"),
    output:
        filtered_barcode_list=f"{output_d}/{{sample}}/filtered_barcode_list.txt",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    resources:
        mem="128G",
        time=24 * 60 - 1,
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_expression')}"
        "python {input.script}"
        " -d {params.d}"
        " --gtf {input.gtf}"
        " --tmp_dir {input.tmp_dir}"
        " -t {threads}"


rule liqa_refgene:
    input:
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        liqa_refgene=f"{output_d}/{{sample}}/liqa.refgene",
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_liqa_refgene')}"
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
        time=ancient(f"{output_d}/time.tsv"),
    output:
        matrix_isoform=f"{output_d}/{{sample}}/matrix_isoform.tsv",
    params:
        d=f"{output_d}/{{sample}}",
    threads: 32
    conda:
        "envs/scnanogps.yaml"
    resources:
        mem="128G",
        time=24 * 60 - 1,
    shell:
        f"{format_time.format_gnu_time_string(process='scNanoGPS_isoforms')}"
        "python {input.script}"
        " --liqa_ref {input.liqa_refgene}"
        " --tmp_dir {input.tmp_dir}"
        " -d {params.d}"
        " -t {threads}"
