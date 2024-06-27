import resource

if len(config) == 0:

    configfile: "Real_config.yaml"


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


use rule * from freddie_smk as freddie_*


use rule * from FLAMES_smk as FLAMES_*


use rule * from scNanoGPS_smk as scNanoGPS_*


output_d = config["outpath"]

tool_names = [
    "freddie",
] + [
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
            f"{output_d}/post/{{sample}}/{{tool_name}}.tsv",
            sample=config["samples"],
            tool_name=tool_names,
        ),
        [
            f"{output_d}/STAR/{sample}/{index_name}/{cell_type}/Aligned.sortedByCoord.out.bam"
            for sample in config["samples"]
            for index_name in ["ENSEMBL", "scFreddie"]
            for cell_type in config["samples"][sample]["rna_srs"]
        ],
    default_target: True


def get_refs(wildcards):
    if wildcards.index_name == "ENSEMBL":
        return dict(
            DNA=config["samples"][wildcards.sample]["DNA"],
            GTF=config["samples"][wildcards.sample]["GTF"],
        )
    elif wildcards.index_name == "scFreddie":
        return dict(
            DNA=config["samples"][wildcards.sample]["DNA"],
            GTF=f"{output_d}/freddie/{{sample}}.isoforms.gtf",
        )


rule STAR_index:
    input:
        fasta=lambda wc: get_refs(wc)["DNA"],
        gtf=lambda wc: get_refs(wc)["GTF"],
    output:
        index=directory(f"{output_d}/STAR/{{sample}}/{{index_name}}_index"),
    threads: 32
    wildcard_constraints:
        index_name="ENSEMBL|scFreddie",
    params:
        tmp_dir=f"{output_d}/STAR/{{sample}}/{{index_name}}_tmp",
    shell:
        "STAR"
        " --runMode genomeGenerate "
        " --genomeDir {output.index} "
        " --genomeFastaFiles {input.fasta} "
        " --sjdbGTFfile {input.gtf}"
        " --runThreadN {threads} "
        " --outTmpDir {params.tmp_dir}"


rule STAR_align:
    input:
        fastq=lambda wc: config["samples"][wc.sample]["rna_srs"][wc.celltype],
        index=f"{output_d}/STAR/{{sample}}/{{index_name}}_index",
    output:
        bam=f"{output_d}/STAR/{{sample}}/{{index_name}}/{{celltype}}/Aligned.sortedByCoord.out.bam",
        SJs=f"{output_d}/STAR/{{sample}}/{{index_name}}/{{celltype}}/SJ.out.tab",
    params:
        prefix=f"{output_d}/STAR/{{sample}}/{{index_name}}/{{celltype}}/",
        limit_n=min(resource.getrlimit(7)) // 32 - 1,
    threads: 32
    shell:
        "STAR"
        " --runMode alignReads"
        " --genomeDir {input.index}"
        " --readFilesIn {input.fastq}"
        " --outFileNamePrefix {params.prefix}"
        " --outSAMtype BAM SortedByCoordinate"
        " --outBAMsortingBinsN {params.limit_n}"
        " --readFilesCommand zcat"
        " --runThreadN {threads}"


rule post_FLAMES:
    input:
        script="py/post_FLAMES.py",
        flames_dir=f"{output_d}/FLAMES/{{sample}}_r{{rate}}",
        cb_tsv=lambda wc: config["samples"][wc.sample]["cb_to_celltype"],
    output:
        tsv=f"{output_d}/post/{{sample}}/FLAMES_r{{rate}}.tsv",
    params:
        gff3=f"{output_d}/FLAMES/{{sample}}_r{{rate}}/isoform_annotated.filtered.gff3",
        csv=f"{output_d}/FLAMES/{{sample}}_r{{rate}}/transcript_count.csv.gz",
    shell:
        "python {input.script}"
        " -gff3 {params.gff3}"
        " -csv {params.csv}"
        " -cb_to_celltypes {input.cb_tsv}"
        " -o {output.tsv}"


rule post_scNanoGPS:
    input:
        script="py/post_scNanoGPS.py",
        gtf=f"{output_d}/scNanoGPS/preprocess/{{sample}}.r{{rate}}.gtf",
        tsv=f"{output_d}/scNanoGPS/{{sample}}_r{{rate}}/matrix_isoform.tsv",
        cb_tsv=lambda wc: config["samples"][wc.sample]["cb_to_celltype"],
    output:
        tsv=f"{output_d}/post/{{sample}}/scNanoGPS_r{{rate}}.tsv",
    shell:
        "python {input.script}"
        " -gtf {input.gtf}"
        " -tsv {input.tsv}"
        " -cb_to_celltypes {input.cb_tsv}"
        " -o {output.tsv}"


rule post_freddie:
    input:
        script="py/post_freddie.py",
        gtf=f"{output_d}/freddie/{{sample}}.isoforms.gtf",
        cb_tsv=lambda wc: config["samples"][wc.sample]["cb_to_celltype"],
    output:
        tsv=f"{output_d}/post/{{sample}}/freddie.tsv",
    shell:
        "python {input.script}"
        " -gtf {input.gtf}"
        " -cb_to_celltypes {input.cb_tsv}"
        " -o {output.tsv}"
