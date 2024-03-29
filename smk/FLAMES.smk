output_d = f"{config['outpath']}/FLAMES"

print("FLAMES:", output_d)


rule isoforms:
    input:
        fastq=f"{output_d}/preprocess/{{sample}}.cb_matched.fastq",
        gtf=lambda wc: config["samples"][wc.sample]["GTF"],
        dna=lambda wc: config["samples"][wc.sample]["DNA"],
        minimap2_tools=f"{output_d}/preprocess/minimap2_tools",
        script=config["FLAMES"]["exec"],
    output:
        directory(
            f"{output_d}/{{sample}}",
        ),
    shell:
        "{input.script}"
        " -a {input.gtf}"
        " -i {input.fastq}"
        " -f {input.dna}"
        " -m {input.minimap2_tools}"
        " -o {output}"


rule match:
    input:
        script=f"{output_d}/match_cell_barcodes",
        fastq=lambda wc: f"{config['samples'][wc.sample]['FASTQ']}",
        cb=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv.gz",
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.cb_stats.tsv",
        fastq=f"{output_d}/preprocess/{{sample}}.cb_matched.fastq",
    params:
        e=2,
        fastq_dir=f"{output_d}/preprocess/{{sample}}_fastqs",
    shell:
        "mkdir -p {params.fastq_dir} && ln -s {input.fastq} {params.fastq_dir}/1.fastq && "
        "{input.script}"
        " {params.fastq_dir}"
        " {output.tsv}"
        " {output.fastq}"
        " {input.cb}"
        " {params.e}"


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
        f"{output_d}/match_cell_barcodes",
    shell:
        """
        g++ -std=c++11 -O2 -I"$(dirname $(dirname  $(which g++)))/include" -o {output} {input}
        """


rule scTagger_extract_bc:
    input:
        tsv=f"{output_d}/preprocess/{{sample}}.lr_bc.tsv.gz",
        cb=lambda wc: config["samples"][wc.sample]["CB"],
    output:
        tsv=f"{output_d}/preprocess/{{sample}}.bc_whitelist.tsv.gz",
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
    shell:
        "scTagger.py extract_lr_bc"
        " -r {input.fastq}"
        " -o {output.tsv}"
        " -t {threads}"


rule minimap2_tools:
    output:
        f"{output_d}/preprocess/minimap2_tools",
    shell:
        """
        mkdir -p {output}
        ln -s $(which minimap2) {output}/minimap2
        ln -s $(which k8) {output}/k8
        ln -s $(which paftools.js) {output}/paftools.js
        """
