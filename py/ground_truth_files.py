#!/usr/bin/env python3
import argparse
import re
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get ground truth table out of MDF files."
    )
    parser.add_argument(
        "-fastq",
        type=str,
        required=True,
        help="Input FASTQ",
    )
    parser.add_argument(
        "-mdfs",
        type=str,
        nargs="+",
        required=True,
        help="Input MDFs",
    )
    parser.add_argument(
        "-o",
        type=str,
        required=True,
        help="Output TSV",
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    rid_to_mid = dict()
    for idx, line in tqdm(enumerate(open(args.fastq)), desc=f"Reading {args.fastq}"):
        if idx % 4 != 0:
            continue
        rid = line.strip().split()[0][1:]
        comments = line.strip().split()[1:]
        for comment in comments:
            if comment.startswith("molecule_id="):
                mid = comment.split("=")[1]
                rid_to_mid[rid] = mid
        assert rid in rid_to_mid
    mid_to_tid = dict()
    mid_to_cb = dict()
    cb_to_cell_lines = dict()
    for mdf in args.mdfs:
        exprmnt = mdf.split("/")[-2]
        for line in tqdm(open(mdf), desc=f"Reading {mdf}"):
            if line[0] != "+":
                continue
            line = line.strip().split("\t")
            mid = line[0][1:]
            cb = "."
            tid = ""
            for comment in line[2].split(";"):
                if comment.startswith("CB"):
                    comment = comment.split("=")
                    if len(comment) == 2:
                        cb = comment[1]
                elif comment.startswith("tid"):
                    tid = comment.split("=")[1]
            mid_to_tid[mid] = tid
            mid_to_cb[mid] = cb
            if cb not in cb_to_cell_lines:
                cb_to_cell_lines[cb] = set()
            cb_to_cell_lines[cb].add(exprmnt)

    outfile = open(args.o, "w+")
    print(
        "\t".join(
            [
                "read_id",
                "molecule_id",
                "transcript_id",
                "cell_barcode",
                "cell_lines",
            ]
        ),
        file=outfile,
    )
    for rid, mid in tqdm(rid_to_mid.items(), desc=f"Writing {args.o}"):
        parent_md = re.split("_|\.", mid)[0]
        record = [
            rid,
            mid,
            mid_to_tid[parent_md],
            mid_to_cb[parent_md],
            ",".join(sorted(cb_to_cell_lines[mid_to_cb[parent_md]])),
        ]
        print("\t".join(record), file=outfile)
    outfile.close()


if __name__ == "__main__":
    main()
