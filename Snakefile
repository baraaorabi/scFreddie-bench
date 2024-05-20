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


rule all:
    input:
        [f"{output_d}/freddie/{s}.isoforms.gtf" for s in config["samples"]],
        [
            f"{output_d}/FLAMES/{s}_r{r}"
            for s in config["samples"]
            for r in config["ref_sample_rates"]
        ],
        [
            f"{output_d}/scNanoGPS/{s}_r{r}/matrix_isoform.tsv"
            for s in config["samples"]
            for r in config["ref_sample_rates"]
        ],
    default_target: True


rule ground_truth_files:
    input:
        fastq=lambda wc: config["samples"][wc.exprmnt]["FASTQ"],
        tsb_mdfs=lambda wc: config["samples"][wc.exprmnt]["tsb_MDFs"],
        fin_mdf=lambda wc: config["samples"][wc.exprmnt]["fin_MDF"],
    output:
        truth_tsv=f"{output_d}/truth/{{exprmnt}}.truth.tsv",
        cb_tsv=f"{output_d}/truth/{{exprmnt}}.cb_to_celltypes.tsv",
    shell:
        "python py/ground_truth_files.py"
        " -fastq {input.fastq}"
        " -tsb_mdfs {input.tsb_mdfs}"
        " -fin_mdf {input.fin_mdf}"
        " -truth_tsv {output.truth_tsv}"
        " -cb_tsv {output.cb_tsv}"
