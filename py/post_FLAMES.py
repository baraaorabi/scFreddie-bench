#!/usr/bin/env python3
import argparse
from collections import defaultdict
import gzip
import os
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post process FLAMES output for accuracy stats."
    )
    parser.add_argument(
        "-gff3",
        type=str,
        required=True,
        help="Input isoform_annotated.filtered.gff3 from FLAMES",
    )
    parser.add_argument(
        "-csv",
        type=str,
        required=True,
        help="Input transcript_count.csv.gz from FLAMES",
    )
    parser.add_argument(
        "-cb_to_celltypes",
        type=str,
        required=True,
        help="Input cb_to_celltypes.tsv file from ground truth",
    )
    parser.add_argument(
        "-o",
        type=str,
        required=True,
        help="Output TSV file",
    )
    args = parser.parse_args()
    return args


def rev_compl(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    args = parse_args()
    print(args)
    cb_to_cts, CTs = get_cb_to_cts(args.cb_to_celltypes)
    tid_to_ct_counts = get_tid_to_ct_counts(cb_to_cts, args.csv)
    tid_to_intervals = get_tid_to_intervals(args.gff3)

    with open(args.o, "w+") as outfile:
        header = ["transcript_id"] + CTs + ["contig", "intervals"]
        print("\t".join(header), file=outfile)
        for tid, (contig, intervals) in tid_to_intervals.items():
            record = [tid]
            record.extend(tid_to_ct_counts[tid][ct] for ct in CTs)
            record.append(contig)
            record.append(",".join(f"{start}-{end}" for start, end in intervals))
            print("\t".join(map(str, record)), file=outfile)


def get_tid_to_intervals(gff3):
    tid_to_intervals: dict[str, tuple[str, list[tuple[int, int]]]] = dict()
    tid = ""
    if not os.path.exists(gff3):
        return tid_to_intervals
    for line in open(gff3):
        if line.startswith("#"):
            continue
        line = line.rstrip("\n").split("\t")
        contig = line[0]
        info = dict(kv.split("=") for kv in line[8].split(";"))
        if line[2] == "transcript":
            tid = info["transcript_id"]
            tid_to_intervals[tid] = (contig, list())
        elif line[2] == "exon":
            assert tid == info["Parent"].split(":")[1]
            start = int(line[3]) - 1
            end = int(line[4])
            tid_to_intervals[tid][1].append((start, end))
    return tid_to_intervals


def get_tid_to_ct_counts(cb_to_cts, csv):
    tid_to_ct_counts = dict()
    if not os.path.exists(csv):
        return tid_to_ct_counts
    with gzip.open(csv, "rt") as f:
        CBs = list()
        for cb in f.readline().rstrip("\n").split(",")[2:]:
            cb_rev = rev_compl(cb)
            if cb <= cb_rev:
                CBs.append(cb)
            else:
                CBs.append(cb_rev)
        for line in tqdm(f, desc=f"Reading {csv}"):
            line = line.rstrip("\n").split(",")
            tid = line[0]
            assert tid not in tid_to_ct_counts
            tid_to_ct_counts[tid] = defaultdict(int)
            for cb, count in zip(CBs, map(int, line[2:])):
                if count == 0:
                    continue
                if cb in cb_to_cts:
                    for ct in cb_to_cts[cb]:
                        tid_to_ct_counts[tid][ct] += count
                else:
                    tid_to_ct_counts[tid]["N"] += count
    return tid_to_ct_counts


def get_cb_to_cts(cb_to_celltypes: str) -> tuple[dict[str, set[str]], list[str]]:
    cb_to_cts: dict[str, set[str]] = dict()
    CT_set = {"N"}
    for line in tqdm(open(cb_to_celltypes), desc=f"Reading {cb_to_celltypes}"):
        cb, cts = line.rstrip("\n").split()
        cts = set(cts.split(","))
        CT_set.update(cts)
        cb_rc = rev_compl(cb)
        if cb <= cb_rc:
            cb_to_cts[cb] = cts
        else:
            cb_to_cts[cb_rc] = cts
    CTs = sorted(CT_set)
    return cb_to_cts, CTs


if __name__ == "__main__":
    main()
