#!/usr/bin/env python3
import argparse
from collections import defaultdict
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post process Freddie output for accuracy stats."
    )
    parser.add_argument(
        "-gtf",
        type=str,
        required=True,
        help="Input annotation GTF from Freddie",
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


def main():
    args = parse_args()

    CTs = get_cts(args.cb_to_celltypes)
    tid_to_ct_counts, tid_to_intervals = process_gtf(gtf=args.gtf)

    with open(args.o, "w+") as outfile:
        header = ["transcript_id"] + CTs + ["contig", "intervals"]
        print("\t".join(header), file=outfile)
        for tid, (contig, intervals) in tid_to_intervals.items():
            record = [tid]
            counts = [tid_to_ct_counts[tid][ct] for ct in CTs]
            if sum(counts) == 0:
                continue
            record.extend(map(str, counts))
            record.append(contig)
            record.append(",".join(f"{start}-{end}" for start, end in intervals))
            print("\t".join(record), file=outfile)


def process_gtf(gtf: str) -> tuple[
    dict[str, dict[str, float]],
    dict[str, tuple[str, list[tuple[int, int]]]],
]:
    tid_to_ct_counts: dict[str, dict[str, float]] = dict()
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
            tid_to_ct_counts[tid] = defaultdict(int)
            for ct in info["cell_types"].split(","):
                if ct == "NA":
                    ct = "N"
                tid_to_ct_counts[tid][ct] = int(info["read_support"])
        elif line[2] == "exon":
            assert tid == info["transcript_id"]
            start = int(line[3]) - 1
            end = int(line[4])
            tid_to_intervals[tid][1].append((start, end))
    return tid_to_ct_counts, tid_to_intervals


def get_cts(cb_to_celltypes: str) -> list[str]:
    CT_set = {"N"}
    for line in tqdm(open(cb_to_celltypes), desc=f"Reading {cb_to_celltypes}"):
        _, cts = line.rstrip("\n").split("\t")
        cts = set(cts.split(","))
        CT_set.update(cts)
    CTs = sorted(CT_set)
    return CTs


if __name__ == "__main__":
    main()
