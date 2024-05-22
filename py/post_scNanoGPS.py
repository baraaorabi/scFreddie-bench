#!/usr/bin/env python3
import argparse
from collections import defaultdict
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post process scNanoGPS output for accuracy stats."
    )
    parser.add_argument(
        "-gtf",
        type=str,
        required=True,
        help="Input annotation GTF used to run scNanoGPS",
    )
    parser.add_argument(
        "-tsv",
        type=str,
        required=True,
        help="Input matrix_isoform.tsv from scNanoGPS",
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

    cb_to_cts, CTs = get_cb_to_cts(args.cb_to_celltypes)
    tid_to_ct_counts = get_tid_to_ct_counts(cb_to_cts=cb_to_cts, tsv=args.tsv)
    tid_to_intervals = get_tid_to_intervals(gtf=args.gtf)

    with open(args.o, "w+") as outfile:
        header = ["transcript_id"] + CTs + ["contig", "intervals"]
        print("\t".join(header), file=outfile)
        for tid, (contig, intervals) in tid_to_intervals.items():
            if tid not in tid_to_ct_counts:
                continue
            record = [tid]
            counts = [round(tid_to_ct_counts[tid][ct]) for ct in CTs]
            if sum(counts) == 0:
                continue
            record.extend(map(str, counts))
            record.append(contig)
            record.append(",".join(f"{start}-{end}" for start, end in intervals))
            print("\t".join(record), file=outfile)


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


def get_tid_to_ct_counts(cb_to_cts, tsv):
    tid_to_ct_counts: dict[str, dict[str, float]] = dict()
    with open(tsv) as f:
        CBs = list()
        for cb in f.readline().rstrip("\n").split("\t")[1:]:
            cb_rev = rev_compl(cb)
            if cb <= cb_rev:
                CBs.append(cb)
            else:
                CBs.append(cb_rev)
        for line in tqdm(f, desc=f"Reading {tsv}"):
            line = line.rstrip("\n").split("\t")
            tid = line[0].split("_")[-1]
            assert len(tid) == 15 and tid.startswith("ENST") and tid[4:].isdigit(), tid
            assert tid not in tid_to_ct_counts
            counts = list(map(float, line[1:]))
            assert len(CBs) == len(counts)
            if sum(counts) == 0:
                continue
            tid_to_ct_counts[tid] = defaultdict(float)
            for cb, count in zip(CBs, counts):
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
        cb, cts = line.rstrip("\n").split("\t")
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
