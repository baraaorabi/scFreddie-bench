import sys
import re

if len(config) == 0:

    configfile: "config.yaml"


module TKSM_smk:
    snakefile:
        "smk/tksm.smk"
    config:
        config["TKSM"]


module freddie_smk:
    snakefile:
        "smk/freddie.smk"
    config:
        config


module FLAMES_smk:
    snakefile:
        "smk/FLAMES.smk"
    config:
        config


module scNanoGPS_smk:
    snakefile:
        "smk/scNanoGPS.smk"
    config:
        config


use rule * from TKSM_smk exclude all as TKSM_*


use rule * from freddie_smk as freddie_*


use rule * from FLAMES_smk as FLAMES_*


use rule * from scNanoGPS_smk as scNanoGPS_*


output_d = config["outpath"]
config["samples"] = dict()
for x in config["TKSM"]["experiments"]:
    final_file = TKSM_smk.exprmnt_final_file(x)
    if final_file.endswith(".fastq"):
        config["samples"][x] = {
            "FASTQ": final_file,
            "DNA": TKSM_smk.get_sample_ref(x, "DNA"),
            "cDNA": TKSM_smk.get_sample_ref(x, "cDNA"),
            "GTF": TKSM_smk.get_sample_ref(x, "GTF"),
            "tsb_MDFs": TKSM_smk.get_source_mdfs(x),
            "fin_MDF": final_file[: -len(".Seq.fastq")] + ".mdf",
            "CB": config["TKSM"]["refs"]["barcodes"]["10x"],
            "truth": f"{output_d}/truth/{x}.truth.tsv",
            "cb_to_celltype": f"{output_d}/truth/{x}.cb_to_celltypes.tsv",
        }

tool_names = ["freddie"] + [
    f"{T}_r{r}"
    for T in [
        "scNanoGPS",
        "FLAMES",
    ]
    for r in config["ref_sample_rates"]
]


rule all:
    input:
        expand(
            f"{output_d}/eval/{{sample}}.{{ext}}",
            sample=config["samples"],
            ext=["pdf", "tsv"],
        ),
    default_target: True


rule ground_truth_files:
    input:
        script="py/ground_truth_files.py",
        fastq=lambda wc: config["samples"][wc.exprmnt]["FASTQ"],
        tsb_mdfs=lambda wc: config["samples"][wc.exprmnt]["tsb_MDFs"],
        fin_mdf=lambda wc: config["samples"][wc.exprmnt]["fin_MDF"],
    output:
        truth_tsv=f"{output_d}/truth/{{exprmnt}}.truth.tsv",
        cb_tsv=f"{output_d}/truth/{{exprmnt}}.cb_to_celltypes.tsv",
    shell:
        "python {input.script}"
        " -fastq {input.fastq}"
        " -tsb_mdfs {input.tsb_mdfs}"
        " -fin_mdf {input.fin_mdf}"
        " -truth_tsv {output.truth_tsv}"
        " -cb_tsv {output.cb_tsv}"


rule post_truth:
    input:
        script="py/post_truth.py",
        bam=f"{output_d}/freddie/preprocess/{{exprmnt}}.sorted.bam",
        truth_tsv=f"{output_d}/truth/{{exprmnt}}.truth.tsv",
        gtf=lambda wc: config["samples"][wc.exprmnt]["GTF"],
    output:
        tsv=f"{output_d}/post/{{exprmnt}}/truth.tsv",
    shell:
        "python {input.script}"
        " -bam {input.bam}"
        " -truth_tsv {input.truth_tsv}"
        " -gtf {input.gtf}"
        " -o {output.tsv}"


rule post_FLAMES:
    input:
        script="py/post_FLAMES.py",
        flames_dir=f"{output_d}/FLAMES/{{exprmnt}}_r{{rate}}",
        cb_tsv=f"{output_d}/truth/{{exprmnt}}.cb_to_celltypes.tsv",
    output:
        tsv=f"{output_d}/post/{{exprmnt}}/FLAMES_r{{rate}}.tsv",
    params:
        gff3=f"{output_d}/FLAMES/{{exprmnt}}_r{{rate}}/isoform_annotated.filtered.gff3",
        csv=f"{output_d}/FLAMES/{{exprmnt}}_r{{rate}}/transcript_count.csv.gz",
    shell:
        "python {input.script}"
        " -gff3 {params.gff3}"
        " -csv {params.csv}"
        " -cb_to_celltypes {input.cb_tsv}"
        " -o {output.tsv}"


rule post_scNanoGPS:
    input:
        script="py/post_scNanoGPS.py",
        gtf=f"{output_d}/scNanoGPS/preprocess/{{exprmnt}}.r{{rate}}.gtf",
        tsv=f"{output_d}/scNanoGPS/{{exprmnt}}_r{{rate}}/matrix_isoform.tsv",
        cb_tsv=f"{output_d}/truth/{{exprmnt}}.cb_to_celltypes.tsv",
    output:
        tsv=f"{output_d}/post/{{exprmnt}}/scNanoGPS_r{{rate}}.tsv",
    shell:
        "python {input.script}"
        " -gtf {input.gtf}"
        " -tsv {input.tsv}"
        " -cb_to_celltypes {input.cb_tsv}"
        " -o {output.tsv}"


rule post_freddie:
    input:
        script="py/post_freddie.py",
        gtf=f"{output_d}/freddie/{{exprmnt}}.isoforms.gtf",
        cb_tsv=f"{output_d}/truth/{{exprmnt}}.cb_to_celltypes.tsv",
    output:
        tsv=f"{output_d}/post/{{exprmnt}}/freddie.tsv",
    shell:
        "python {input.script}"
        " -gtf {input.gtf}"
        " -cb_to_celltypes {input.cb_tsv}"
        " -o {output.tsv}"


rule evaluate:
    input:
        script="py/evaluate.py",
        tsvs=[f"{output_d}/post/{{sample}}/{t}.tsv" for t in tool_names],
        truth=f"{output_d}/post/{{sample}}/truth.tsv",
    output:
        pdf=f"{output_d}/eval/{{sample}}.pdf",
        tsv=f"{output_d}/eval/{{sample}}.tsv",
    params:
        names=tool_names,
    shell:
        "python {input.script}"
        " -sample {wildcards.sample}"
        " -tsvs {input.tsvs}"
        " -truth {input.truth}"
        " -names {params.names}"
        " -pdf {output.pdf}"
        " -tsv {output.tsv}"
