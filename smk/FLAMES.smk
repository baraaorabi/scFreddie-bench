import sys
import re
from collections import Counter


preproc_d = f"{config['outpath']}/FLAMES"
FLAME_path = config["FLAME_path"]


rule make_match_cell_barcodes:
    input:
        expand(
            f"{FLAME_path}/src/{f}",
            f=[
                "ssw/ssw_cpp.cpp",
                "ssw/ssw.c",
                "match_cell_barcode.cpp",
                "kseq.h",
                "edit_dist.cp",
            ],
        ),
    output:
        f"{FLAME_path}/bin/match_cell_barcodes",
    shell:
        """
        g++ -std=c++11 -O2 -I"$(dirname $(dirname  $(which g++)))/include" -o {output} {input}
        """

rule match_cell_barcodes:
    input:
        script=f"{FLAME_path}/bin/match_cell_barcodes",
        fastq=lambda wc: f"{config['samples'][wc.sample]['fastq']}",
        barcode=lambda wc: f"{config['samples'][wc.sample]['barcode']}",
    output:
        tsv=f"{preproc_d}/{sample}.cb_stats.tsv",
        fastq=f"{preproc_d}/{sample}.cb_matched.fastq",
    params:
        e=2,
        fastq_dir=f"{preproc_d}/{sample}_fastqs",
    shell:
        "mkdir -p {params.fastq_dir} && ln -s {input.fastq} {params.fastq_dir}/1.fastq"
        " && {input.script} {params.fastq_dir} {output.tsv} {output.fastq} {input.barcode} {params.e}"
