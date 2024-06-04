#!/usr/bin/env python3
import argparse
from glob import glob

import cgranges
import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt


count_range = [3, 5, 9] + list(range(15, 51, 5))


def parse_args():
    parser = argparse.ArgumentParser(description="Plot tool performance.")
    parser.add_argument(
        "-sample",
        type=str,
        required=True,
        help="Sample name",
    )
    parser.add_argument(
        "-tsvs",
        type=str,
        required=True,
        nargs="+",
        help="Input TSVs for the tools",
    )
    parser.add_argument(
        "-truth",
        type=str,
        required=True,
        help="Input TSV for the truth",
    )
    parser.add_argument(
        "-names",
        type=str,
        required=True,
        nargs="+",
        help="Tool names",
    )
    parser.add_argument(
        "-out",
        type=str,
        required=True,
        help="Output filename prefix. Will generate <output>.tsv, <output>.pdf, <output>.png, and <output>.svg",
    )

    args = parser.parse_args()
    return args


def parse_intervals(intervals):
    return [tuple(map(int, interval.split("-"))) for interval in intervals.split(",")]


def parse_tsv(tsv):
    df = pd.read_csv(tsv, sep="\t")
    CTs = [s for s in df.columns[1:-2]]
    df["contig"] = df["contig"].astype(str)
    df["intervals"] = df["intervals"].apply(parse_intervals)
    df = df[df.apply(lambda x: len(x["intervals"]) > 1, axis=1)]
    df = pd.melt(
        df,
        id_vars=["transcript_id", "contig", "intervals"],
        value_vars=CTs,
        var_name="cell_type",
        value_name="count",
    )
    df = df[
        [
            "transcript_id",
            "cell_type",
            "count",
            "contig",
            "intervals",
        ]
    ]
    df = df[df["count"] > 0]
    df = df.sort_values("count", ascending=False)
    df = df.reset_index(drop=True)
    return df


def get_overlaps(pred_record):
    tidxs = set()
    contig = pred_record["contig"]
    cell_type = pred_record["cell_type"]
    for start, end in pred_record["intervals"]:
        for _, _, tidx in truth_cg.overlap(f"{cell_type}-{contig}", start, end):
            tidxs.add(tidx)
    return sorted(tidxs)


def fuzzy_match(A, B, slack=20):
    if len(A) != len(B):
        return False
    for (A_s, A_e), (B_s, B_e) in zip(A, B):
        if abs(B_s - A_s) > slack or abs(B_e - A_e) > slack:
            return False
    return True


def get_intronic_matches(pred_record):
    tidxs = list()
    for tidx in get_overlaps(pred_record):
        assert truth_df["contig"][tidx] == pred_record["contig"]
        A = truth_df["intervals"][tidx]
        B = pred_record["intervals"]
        A_introns = [(e1, s2) for (_, e1), (s2, _) in zip(A[:-1], A[1:])]
        B_introns = [(e1, s2) for (_, e1), (s2, _) in zip(B[:-1], B[1:])]
        if fuzzy_match(A_introns, B_introns):
            tidxs.append(tidx)
    return tidxs


def set_truth_df(tsv):
    global truth_df, truth_cg
    truth_df = parse_tsv(tsv)

    truth_cg = cgranges.cgranges()
    for idx, (contig, intervals, cell_type) in enumerate(
        zip(
            truth_df["contig"],
            truth_df["intervals"],
            truth_df["cell_type"],
        )
    ):
        for start, end in intervals:
            truth_cg.add(f"{cell_type}-{contig}", start, end, idx)
    truth_cg.index()


def get_pred_dfs(tsvs, names):
    pred_dfs = dict()
    for tool, tsv in tqdm(zip(names, tsvs), total=len(tsvs), desc="Parsing TSVs"):
        if tool == "truth":
            continue
        try:
            pred_df = parse_tsv(tsv)
            pred_df["hits"] = pred_df.apply(get_intronic_matches, axis=1)
        except:
            continue
        pred_dfs[tool] = pred_df
    return pred_dfs


def plot(curves, outprefix, sample):
    fig, axes = plt.subplots(1, 3, figsize=(10, 3), sharey=True)
    txt = fig.suptitle(f"Evaluation for {sample}")
    axes[0].set_xlabel("Min truth count")
    axes[0].set_ylabel("Recall")
    axes[1].set_xlabel("Min tool support")
    axes[1].set_ylabel("Precision")
    axes[2].set_xlabel("Min truth count and tool support")
    axes[2].set_ylabel("F1 score")

    for tool, (R_curve, P_curve, f1_curve) in curves.items():
        kw = styling(tool)
        axes[0].plot(count_range[: len(R_curve)], R_curve, label=tool, **kw)
        axes[1].plot(count_range[: len(P_curve)], P_curve, label=tool, **kw)
        axes[2].plot(count_range[: len(f1_curve)], f1_curve, label=tool, **kw)
    for ax in axes:
        ax.set_xticks(count_range)
        ax.set_xticklabels(count_range)
        ax.set_ylim(0, 1)
    fig.tight_layout()
    lgd = plt.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))
    for ext in ["pdf", "png", "svg"]:
        plt.savefig(
            f"{outprefix}.{ext}",
            bbox_extra_artists=(lgd, txt),
            bbox_inches="tight",
        )


def styling(tool):
    D = dict()
    if tool == "scFreddie":
        D["linestyle"] = "-."
        D["marker"] = "s"
        D["color"] = "#1b9e77"
        D["zorder"] = 3
    elif tool.startswith("scNanoGPS"):
        D["linestyle"] = "-"
        D["marker"] = "x"
        D["color"] = "#d95f02"
        D["zorder"] = 2
        r = float(tool.split("_")[1][1:])
        D["alpha"] = 0.4 + 0.6 * (r - 0.01) / 0.99
    elif tool.startswith("FLAMES"):
        D["linestyle"] = "--"
        D["marker"] = "o"
        D["color"] = "#7570b3"
        D["zorder"] = 1
        r = float(tool.split("_")[1][1:])
        D["alpha"] = 0.4 + 0.6 * (r - 0.01) / 0.99
    D["markersize"] = 3
    return D


def get_curves(pred_dfs):
    curves = dict()
    for tool, pred_df in pred_dfs.items():
        R = np.zeros(len(truth_df), dtype=bool)
        for tidxs in pred_df["hits"]:
            for tidx in tidxs:
                R[tidx] = True
        R_curve = list()
        for count in count_range:
            num = np.sum(R[truth_df["count"] >= count])
            den = np.sum(truth_df["count"] >= count)
            if den == 0:
                break
            R_curve.append(num / den)
        P = np.array(pred_df["hits"].apply(lambda x: len(x) > 0), dtype=bool)
        P_curve = list()
        for count in count_range:
            num = np.sum(P[pred_df["count"] >= count])
            den = np.sum(pred_df["count"] >= count)
            if den == 0:
                break
            P_curve.append(num / den)
        l = min(len(R_curve), len(P_curve))
        f1_curve = (
            2
            * (np.array(P_curve[:l]) * np.array(R_curve[:l]))
            / (np.array(P_curve[:l]) + np.array(R_curve[:l]))
        )
        curves[tool] = (R_curve, P_curve, f1_curve)
    return curves


def output_curves(curves, outprefix):
    tsv = outprefix + ".tsv"
    with open(tsv, "w") as f:
        f.write("tool\tmin_count\trecall\tprecision\tf1\n")
        for tool, (R_curve, P_curve, f1_curve) in curves.items():
            for idx, min_count in enumerate(count_range):
                r = "NA" if idx >= len(R_curve) else f"{R_curve[idx]:.2f}"
                p = "NA" if idx >= len(P_curve) else f"{P_curve[idx]:.2f}"
                f1 = "NA" if idx >= len(f1_curve) else f"{f1_curve[idx]:.2f}"
                f.write(f"{tool}\t{min_count}\t{r}\t{p}\t{f1}\n")
        f.close()


def main():
    args = parse_args()

    set_truth_df(args.truth)
    pred_dfs = get_pred_dfs(args.tsvs, args.names)
    curves = get_curves(pred_dfs)
    output_curves(curves, args.out)
    plot(curves, args.out, args.sample)


if __name__ == "__main__":
    main()
