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


output_d = config["outpath"]
tksm_d = f"{output_d}/TS"
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
        [
            path.format(e=e)
            for e in fastq_exprmnts
            for path in [
                f"{tksm_d}/{{e}}.fastq",
                f"{tksm_d}/{{e}}.tsv",
                f"{output_d}/{{e}}/{{e}}.sorted.bam",
                f"{output_d}/scTagger/{{e}}/{{e}}.bc_whitelist.tsv.gz",
                f"{output_d}/rname_to_celltypes/{{e}}.tsv",
            ]
        ],
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
        "cp {input.fastq} {output.fastq} && "
        "python py/ground_truth_files.py "
        " -fastq {output.fastq}"
        " -mdfs {input.mdfs}"
        " -o {output.tsv}"


rule main_minimap2:
    input:
        reads=f"{tksm_d}/{{sample}}.fastq",
        genome=lambda wc: TS_smk.get_sample_ref(wc.sample, "DNA"),
    output:
        bam=f"{output_d}/{{sample}}/{{sample}}.sorted.bam",
        bai=f"{output_d}/{{sample}}/{{sample}}.sorted.bam.bai",
    threads: 32
    resources:
        mem="128G",
        time=24 * 60 - 1,
    shell:
        "minimap2 -a -x splice -t {threads} {input.genome} {input.reads} | "
        "  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam - > {output.bam} && "
        "  samtools index {output.bam}"


rule main_scTagger_match:
    input:
        lr_tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
        wl_tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
    output:
        lr_tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
    threads: 32
    shell:
        "scTagger.py match_trie"
        " -lr {input.lr_tsv}"
        " -sr {input.wl_tsv}"
        " -o {output.lr_tsv}"
        " -t {threads}"


rule main_scTagger_extract_bc:
    input:
        tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
        wl=config["refs"]["barcodes"]["10x"],
    output:
        tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.bc_whitelist.tsv.gz",
    conda:
        "Snakemake-envs/sctagger.yaml"
    shell:
        "scTagger.py extract_sr_bc_from_lr"
        " -i {input.tsv}"
        " -wl {input.wl}"
        " -o {output.tsv}"


rule main_scTagger_lr_seg:
    input:
        reads=f"{tksm_d}/{{sample}}.fastq",
    output:
        tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.lr_bc.tsv.gz",
    threads: 32
    shell:
        "scTagger.py extract_lr_bc"
        " -r {input.reads}"
        " -o {output.tsv}"
        " -t {threads}"


rule rname_to_celltypes:
    input:
        lr_matches_tsv=f"{output_d}/scTagger/{{sample}}/{{sample}}.lr_matches.tsv.gz",
        truth_tsv=f"{tksm_d}/{{sample}}.tsv",
    output:
        tsv=f"{output_d}/rname_to_celltypes/{{sample}}.tsv",
    shell:
        "python py/rname_to_celltypes.py"
        " -lr {input.lr_matches_tsv}"
        " -truth {input.truth_tsv}"
        " -o {output.tsv}"
