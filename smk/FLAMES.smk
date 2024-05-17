import gzip
import format_time

output_d = f"{config['outpath']}/FLAMES"

minimap2_tools_list = [
    "minimap2",
    "k8",
    "paftools.js",
]


rule make_time:
    output:
        time=f"{output_d}/time.tsv",
    run:
        with open(output.time, "w+") as fout:
            print(format_time.header_str, file=fout)


rule isoforms:
    input:
        fastq=f"{output_d}/preprocess/{{sample}}.cb_matched.fastq.gz",
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
        dna=lambda wc: config["samples"][wc.sample]["DNA"],
        minimap2_tools=[f"{output_d}/exec/{t}" for t in minimap2_tools_list],
        script=config["FLAMES"]["exec"],
        time=ancient(f"{output_d}/time.tsv"),
    params:
        minimap2_tools=f"{output_d}/exec",
    output:
        directory(
            f"{output_d}/{{sample}}",
        ),
    conda:
        "envs/flames.yaml"
    resources:
        mem="128G",
        time=72 * 60 - 1,
    threads: 12
    shell:
        f"{format_time.format_gnu_time_string(process='FLAMES_isoforms')}"
        "{input.script}"
        " -a {input.gtf}"
        " -i {input.fastq}"
        " -f {input.dna}"
        " -m {params.minimap2_tools}"
        " -o {output}"


rule match:
    input:
        script=f"{output_d}/exec/match_cell_barcodes",
        cb=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv",
        fastq=f"{output_d}/preprocess/{{sample}}_fastqs/1.fastq",
        time=ancient(f"{output_d}/time.tsv"),
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.cb_stats.tsv",
        fastq=f"{output_d}/preprocess/{{sample}}.cb_matched.fastq.gz",
    params:
        e=2,
        fastq_dir=f"{output_d}/preprocess/{{sample}}_fastqs",
    conda:
        "envs/flames.yaml"
    resources:
        mem="128G",
        time=24 * 60 - 1,
    shell:
        f"{format_time.format_gnu_time_string(process='FLAMES_match')}"
        "{input.script}"
        " {params.fastq_dir}"
        " {output.tsv}"
        " {output.fastq}"
        " {input.cb}"
        " {params.e}"


rule ln_fastq:
    input:
        fastq=lambda wc: f"{config['samples'][wc.sample]['FASTQ']}",
    output:
        fastq=f"{output_d}/preprocess/{{sample}}_fastqs/1.fastq",
    shell:
        "ln -s $(readlink -f {input.fastq}) {output.fastq}"


rule make_match:
    input:
        [
            f"{config['FLAMES']['src']}/{f}"
            for f in [
                "ssw/ssw_cpp.cpp",
                "ssw/ssw.c",
                "match_cell_barcode.cpp",
                "kseq.h",
                "edit_dist.cpp",
            ]
        ],
    output:
        f"{output_d}/exec/match_cell_barcodes",
    conda:
        "envs/flames.yaml"
    shell:
        """
        g++ -std=c++11 -O2 -lz -I"$(dirname $(dirname  $(which g++)))/include" -o {output} {input}
        """


rule unzip_bc:
    input:
        tsv=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv.gz",
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv",
    run:
        with gzip.open(input.tsv, "rt") as f_in, open(output.tsv, "w+") as f_out:
            for line in f_in:
                bc = line.split("\t")[0]
                f_out.write(f"{bc}-1\n")


rule scTagger_extract_bc:
    input:
        tsv=f"{output_d}/preprocess/{{sample}}.lr_bc.tsv.gz",
        cb=lambda wc: config["samples"][wc.sample]["CB"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv.gz",
    conda:
        "envs/sctagger.yaml"
    shell:
        f"{format_time.format_gnu_time_string(process='FLAMES_scTagger_extract_bc')}"
        "scTagger.py extract_sr_bc_from_lr"
        " -i {input.tsv}"
        " -wl {input.cb}"
        " -o {output.tsv}"


rule scTagger_lr_seg:
    input:
        fastq=lambda wc: config["samples"][wc.sample]["FASTQ"],
        time=ancient(f"{output_d}/time.tsv"),
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.lr_bc.tsv.gz",
    threads: 32
    conda:
        "envs/sctagger.yaml"
    shell:
        f"{format_time.format_gnu_time_string(process='FLAMES_scTagger_lr_seg')}"
        "scTagger.py extract_lr_bc"
        " -r {input.fastq}"
        " -o {output.tsv}"
        " -t {threads}"


rule minimap2_tools:
    output:
        [f"{output_d}/exec/{t}" for t in minimap2_tools_list],
    params:
        out_d=f"{output_d}/exec",
    conda:
        "envs/flames.yaml"
    shell:
        """
        ln -s $(which minimap2) {params.out_d}/minimap2
        ln -s $(which k8) {params.out_d}/k8
        ln -s $(which paftools.js) {params.out_d}/paftools.js
        """
