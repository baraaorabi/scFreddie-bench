#!/usr/bin/env python3
import argparse
from collections import defaultdict
import gzip
import re
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get cell typing of reads using ground truth and scTagger matching."
    )
    parser.add_argument(
        "-lr_matches_tsv",
        type=str,
        required=True,
        help="Input: scTagger matches TSV",
    )
    parser.add_argument(
        "-truth_tsv",
        type=str,
        required=True,
        help="Input: truth TSV",
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
    cb_to_celltypes = defaultdict(set)
    with open(args.truth_tsv, "r") as f:
        f.readline()
        for line in tqdm(f, desc=f"Reading {args.truth_tsv}"):
            rname, _, _, cb, ct = line.strip().split()
            if cb == ".":
                continue
            cb_to_celltypes[cb].update(ct.split(","))

    rname_to_celltypes = defaultdict(set)
    with gzip.open(args.lr_matches_tsv, "rt") as f:
        for line in tqdm(f, desc=f"Reading {args.lr_matches_tsv}"):
            rname, _, _, _, barcodes = line.strip().split()
            for bc in barcodes.split(","):
                rname_to_celltypes[rname].update(cb_to_celltypes[bc])

    with open(args.o, "w+") as f:
        for rname, celltypes in tqdm(rname_to_celltypes.items(), desc=f"Writing {args.o}"):
            f.write(f"{rname}\t{','.join(celltypes)}\n")


if __name__ == "__main__":
    main()
