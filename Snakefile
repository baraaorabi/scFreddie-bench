import sys
import re

if len(config) == 0:

    configfile: "config.yaml"


module TS_smk:
    snakefile:
        "smk/tksm.smk"
    config:
        config


### Import TKSM Snakemake rules
use rule * from TS_smk exclude all


outpath = config["outpath"]
preproc_d = f"{outpath}/preprocess"
tksm_d = f"{outpath}/TS"
exprmnts_re = "|".join([re.escape(x) for x in config["TS_experiments"]])
fastq_exprmnts = [
    x
    for x in config["TS_experiments"]
    if TS_smk.exprmnt_final_file(x).endswith(".fastq")
]
fastq_exprmnts_re = "|".join([re.escape(x) for x in fastq_exprmnts])


def get_source_mdfs(exprmnt):
    first_step = config["TS_experiments"][exprmnt]["pipeline"][0]
    rule_name = list(first_step)[0]
    first_step = first_step[rule_name]
    if rule_name == "Mrg":
        return sorted(
            {mdf for source in first_step["sources"] for mdf in get_source_mdfs(source)}
        )
    if rule_name == "Tsb":
        return [f"{tksm_d}/{exprmnt}/Tsb.mdf"]


rule all:
    input:
        expand(
            f"{tksm_d}/{{exprmnt}}.{{ext}}",
            exprmnt=fastq_exprmnts,
            ext=["fastq", "tsv"],
        ),
    default_target: True


rule ground_truth_files:
    input:
        fastq=lambda wc: TS_smk.exprmnt_final_file(wc.exprmnt),
        mdfs=lambda wc: get_source_mdfs(wc.exprmnt),
    output:
        fastq=f"{tksm_d}/{{exprmnt}}.fastq",
        tsv=f"{tksm_d}/{{exprmnt}}.tsv",
    wildcard_constraints:
        exprmnt=fastq_exprmnts_re,
    shell:
        "cp {input.fastq} {output.fastq}"
        " && python py/ground_truth_files.py -fastq {output.fastq} -mdfs {input.mdfs} -o {output.tsv}"
