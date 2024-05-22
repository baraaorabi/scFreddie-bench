#!/usr/bin/env python3
import typing

import pysam
import pandas as pd
from tqdm import tqdm
import cgranges
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get detectable list of transcripts based on their mapping."
    )
    parser.add_argument(
        "-bam",
        type=str,
        required=True,
        help="BAM file of TKSM reads",
    )
    parser.add_argument(
        "-truth_tsv",
        type=str,
        required=True,
        help="Ground truth TSV file of TKSM reads",
    )
    parser.add_argument(
        "-gtf",
        type=str,
        required=True,
        help="Annotation GTF",
    )
    parser.add_argument(
        "-o",
        type=str,
        required=True,
        help="Output TSV",
    )

    args = parser.parse_args()
    return args


def get_gtf_data(gtf: str):
    tid_to_gidx: dict[str, int] = dict()
    genes: dict[typing.Union[int, str], tuple[int, str]] = dict()
    gene_tree = cgranges.cgranges()
    for line in tqdm(open(gtf), desc=f"Reading {gtf}"):
        if line[0] == "#":
            continue
        line = line.strip("\n").split("\t")
        if not line[2] in ["gene", "transcript"]:
            continue
        info = line[8]
        info = [x.strip().split(" ") for x in info.strip(";").split(";")]
        info = {x[0]: x[1].strip('"') for x in info}
        gid = info["gene_id"]
        if not gid in genes:
            gidx = len(genes) // 2
            genes[gidx] = (gidx, gid)
            genes[gid] = (gidx, gid)
        gidx = genes[gid][0]
        if line[2] == "gene":
            c, s, e = line[0], int(line[3]), int(line[4])
            gene_tree.add(c, s, e, genes[gid][0])
        elif line[2] == "transcript":
            tid_to_gidx[info["transcript_id"]] = gidx
    gene_tree.index()
    return tid_to_gidx, gene_tree


def add_mapping_info(bam, tid_to_gidx, gene_tree, reads):
    rid_to_ridx = {rid: idx for idx, rid in enumerate(reads.read_id)}
    maps_to_gene = [False] * len(rid_to_ridx)
    for aln in tqdm(pysam.AlignmentFile(bam), desc=f"Reading {bam}"):
        if aln.is_supplementary or aln.is_secondary or aln.is_unmapped:
            continue
        ridx = rid_to_ridx[aln.query_name]
        tid = reads.transcript_id[ridx]
        gidx = tid_to_gidx[tid]
        c = aln.reference_name
        s = aln.reference_start
        e = aln.reference_end
        for _, _, i in gene_tree.overlap(c, s, e):
            if i == gidx:
                maps_to_gene[ridx] = True
                break
    reads["maps_to_gene"] = maps_to_gene


def get_transcripts(reads, tid_to_intervals):
    CTs = sorted(set(reads["cell_lines"]))
    tids = sorted(set(reads["transcript_id"]))
    transcripts = pd.DataFrame(
        columns=[
            "transcript_id",
        ]
        + CTs
        + ["contig", "intervals"],
    )
    transcripts["transcript_id"] = tids
    for ct in CTs:
        transcripts[ct] = 0
    transcripts["contig"] = ""
    transcripts["intervals"] = ""
    transcripts = transcripts.set_index("transcript_id", verify_integrity=True)

    for tid, (contig, intervals) in tid_to_intervals.items():
        if not tid in transcripts.index:
            continue
        transcripts.at[tid, "contig"] = contig
        transcripts.at[tid, "intervals"] = ",".join([f"{s}-{e}" for s, e in intervals])

    for _, read in tqdm(reads.iterrows(), total=len(reads), desc="Counting reads"):
        if read["maps_to_gene"]:
            ct = read["cell_lines"]
            tid = read["transcript_id"]
            transcripts.at[tid, ct] += 1
    transcripts = transcripts.sort_values(["contig", "intervals"]).reset_index()
    return transcripts


def get_tid_to_intervals(gtf):
    tid_to_intervals: dict[str, tuple[str, list[tuple[int, int]]]] = dict()
    tid = ""
    for line in open(gtf):
        if line.startswith("#"):
            continue
        line = line.rstrip("\n").split("\t")
        contig = line[0]
        info = line[8]
        info = [x.strip().split(" ") for x in info.strip(";").split(";")]
        info = {x[0]: x[1].strip('"') for x in info}
        if line[2] == "transcript":
            tid = info["transcript_id"]
            tid_to_intervals[tid] = (contig, list())
        elif line[2] == "exon":
            assert tid == info["transcript_id"]
            start = int(line[3]) - 1
            end = int(line[4])
            tid_to_intervals[tid][1].append((start, end))
    return tid_to_intervals


def main():
    args = parse_args()
    tid_to_gidx, gene_tree = get_gtf_data(gtf=args.gtf)
    reads = (
        pd.read_csv(
            args.truth_tsv,
            sep="\t",
        )
        .drop(
            ["cell_barcode"],
            axis=1,
        )
        .fillna("N")
    )

    add_mapping_info(
        bam=args.bam,
        tid_to_gidx=tid_to_gidx,
        gene_tree=gene_tree,
        reads=reads,
    )

    reads["cell_lines"] = reads["cell_lines"].str.split(",")
    reads = reads.explode("cell_lines", ignore_index=True)
    transcripts = get_transcripts(
        reads=reads,
        tid_to_intervals=get_tid_to_intervals(gtf=args.gtf),
    )
    transcripts.to_csv(args.o, sep="\t", index=False)


if __name__ == "__main__":
    main()
