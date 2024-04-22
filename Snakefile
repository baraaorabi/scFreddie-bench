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


use rule * from TKSM_smk exclude all as TKSM_*


use rule * from freddie_smk as freddie_*


use rule * from FLAMES_smk as FLAMES_*


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
            "GFF3": TKSM_smk.get_sample_ref(x, "GFF3"),
            "MDFs": TKSM_smk.get_source_mdfs(x),
            "CB": config["TKSM"]["refs"]["barcodes"]["10x"],
            "truth": f"{output_d}/truth/{x}.tsv",
        }


rule all:
    input:
        [f for s in config["samples"].values() for f in s.values()],
        [f"{output_d}/freddie/{s}.isoforms.gtf" for s in config["samples"]],
        [f"{output_d}/FLAMES/{s}" for s in config["samples"]],
    default_target: True


rule ground_truth_files:
    input:
        fastq=lambda wc: config["samples"][wc.exprmnt]["FASTQ"],
        mdfs=lambda wc: config["samples"][wc.exprmnt]["MDFs"],
    output:
        truth=f"{output_d}/truth/{{exprmnt}}.tsv",
    shell:
        "python py/ground_truth_files.py"
        " -fastq {input.fastq}"
        " -mdfs {input.mdfs}"
        " -o {output.truth}"
