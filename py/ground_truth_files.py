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
        "-tsb_mdfs",
        type=str,
        nargs="+",
        required=True,
        help="Input Tsb MDFs",
    )
    parser.add_argument(
        "-fin_mdf",
        type=str,
        required=True,
        help="Input final MDF before sequencing",
    )
    parser.add_argument(
        "-truth_tsv",
        type=str,
        required=True,
        help="Output truth TSV",
    )
    parser.add_argument(
        "-cb_tsv",
        type=str,
        required=True,
        help="Output cell barcode to cell line/type TSV",
    )

    args = parser.parse_args()
    return args


def process_Tsb_mdfs(
    tsb_mdfs: list[str],
):
    tmid_to_tid: dict[str, str] = dict()
    tmid_to_cb: dict[str, str] = dict()
    cb_to_cell_lines: dict[str, set[str]] = {".": set()}
    for mdf in tsb_mdfs:
        for line in tqdm(open(mdf), desc=f"Reading {mdf}"):
            if line[0] != "+":
                continue
            line = line.rstrip("\n").split("\t")
            tmid = line[0][1:]
            cell_line = tmid.split(":")[0]
            cb = "."
            tid = ""
            for comment in line[2].split(";"):
                if comment.startswith("CB"):
                    comment = comment.split("=")
                    if len(comment) == 2:
                        cb = comment[1]
                elif comment.startswith("tid"):
                    tid = comment.split("=")[1]
            tmid_to_tid[tmid] = tid
            tmid_to_cb[tmid] = cb
            if cb not in cb_to_cell_lines:
                cb_to_cell_lines[cb] = set()
            cb_to_cell_lines[cb].add(cell_line)
    return tmid_to_tid, tmid_to_cb, cb_to_cell_lines


def process_fin_mdf(fin_mdf: str, tmid_to_cb: dict[str, str]):
    fmid_to_cb: dict[str, str] = dict()
    fmid_to_tmid: dict[str, str] = dict()
    for line in tqdm(open(fin_mdf), desc=f"Reading {fin_mdf}"):
        line = line.rstrip("\n").split("\t")
        if line[0][0] == "+":
            fmid = line[0][1:]
            tmid = re.split(r"[_.]", fmid)[0]
            tcb = tmid_to_cb[tmid]
            fmid_to_cb[fmid] = "."
            fmid_to_tmid[fmid] = tmid
            continue
        contig, start, end, _, _ = line
        if contig == tcb and (len(tcb) - (int(end) - int(start)) <= 2):
            fmid_to_cb[fmid] = tcb
    return fmid_to_cb, fmid_to_tmid


def process_fastq(
    fastq: str,
    fmid_to_tmid: dict[str, str],
    fmid_to_cb: dict[str, str],
):
    rid_to_fmid: dict[str, str] = dict()
    for idx, line in tqdm(enumerate(open(fastq)), desc=f"Reading {fastq}"):
        if idx % 4 != 0:
            continue
        rid = line.strip().split()[0][1:]
        comments = line.strip().split()[1:]
        for comment in comments:
            if comment.startswith("molecule_id="):
                fmid = comment.split("=")[1]
                assert fmid in fmid_to_tmid and fmid in fmid_to_cb
                rid_to_fmid[rid] = fmid
        assert rid in rid_to_fmid
    return rid_to_fmid


def write_truth_tsv(
    truth_tsv: str,
    rid_to_fmid: dict[str, str],
    fmid_to_tmid: dict[str, str],
    fmid_to_cb: dict[str, str],
    tmid_to_tid: dict[str, str],
    cb_to_cell_lines: dict[str, set[str]],
):
    outfile = open(truth_tsv, "w+")
    print(
        "\t".join(
            [
                "read_id",
                "transcript_id",
                "cell_barcode",
                "cell_lines",
            ]
        ),
        file=outfile,
    )
    for rid, fmid in tqdm(rid_to_fmid.items(), desc=f"Writing {truth_tsv}"):
        tmid = fmid_to_tmid[fmid]
        cb = fmid_to_cb[fmid]
        cell_lines = cb_to_cell_lines[cb]
        record = [
            rid,
            tmid_to_tid[tmid],
            cb,
            ",".join(sorted(cell_lines)),
        ]
        print("\t".join(record), file=outfile)
    outfile.close()


def write_cb_tsv(
    cb_tsv: str,
    cb_to_cell_lines: dict[str, set[str]],
):
    outfile = open(cb_tsv, "w+")
    for cb, cell_lines in cb_to_cell_lines.items():
        print(
            cb,
            ",".join(sorted(cell_lines)),
            file=outfile,
        )
    outfile.close()


def main():
    args = parse_args()

    tmid_to_tid, tmid_to_cb, cb_to_cell_lines = process_Tsb_mdfs(
        tsb_mdfs=args.tsb_mdfs,
    )
    fmid_to_cb, fmid_to_tmid = process_fin_mdf(
        fin_mdf=args.fin_mdf,
        tmid_to_cb=tmid_to_cb,
    )
    rid_to_fmid = process_fastq(
        fastq=args.fastq,
        fmid_to_tmid=fmid_to_tmid,
        fmid_to_cb=fmid_to_cb,
    )

    write_truth_tsv(
        truth_tsv=args.truth_tsv,
        rid_to_fmid=rid_to_fmid,
        fmid_to_tmid=fmid_to_tmid,
        fmid_to_cb=fmid_to_cb,
        tmid_to_tid=tmid_to_tid,
        cb_to_cell_lines=cb_to_cell_lines,
    )
    write_cb_tsv(
        cb_tsv=args.cb_tsv,
        cb_to_cell_lines=cb_to_cell_lines,
    )


if __name__ == "__main__":
    main()
