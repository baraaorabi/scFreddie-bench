#!/usr/bin/env python3
from dataclasses import dataclass, field
import multiprocessing
from collections import Counter
import enum
from functools import total_ordering, partial
from itertools import groupby
from typing import Generator
import argparse

import pandas as pd
import pysam
from tqdm import tqdm


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


@dataclass
class Intervals:
    contig: str = ""
    intervals: list[tuple[int, int]] = field(default_factory=list)


def get_gtf_data(gtf: str) -> dict[str, Intervals]:
    transcript_intervals: dict[str, Intervals] = dict()
    for line in tqdm(open(gtf), desc=f"Reading {gtf}"):
        if line[0] == "#":
            continue
        line = line.strip("\n").split("\t")
        if not line[2] in ["exon", "transcript"]:
            continue
        info = line[8]
        info = [x.strip().split(" ") for x in info.strip(";").split(";")]
        info = {x[0]: x[1].strip('"') for x in info}
        tid = info["transcript_id"]
        contig, start, end = line[0], int(line[3]) - 1, int(line[4]) - 1
        match line[2]:
            case "transcript":
                assert tid not in transcript_intervals
                transcript_intervals[tid] = Intervals(
                    contig=contig,
                )
            case "exon":
                transcript_intervals[tid].intervals.append((start, end))
    transcript_intervals = {
        k: v for k, v in transcript_intervals.items() if len(v.intervals) > 1
    }
    return transcript_intervals


def generate_introns(
    intervals: list[tuple[int, int]]
) -> Generator[tuple[int, int], None, None]:
    for (_, e1), (s2, _) in zip(intervals[:-1], intervals[1:]):
        yield (e1, s2)


def close_enough(t_intron, r_intron, slack=10):
    t_start, t_end = t_intron
    r_start, r_end = r_intron
    if abs(t_start - r_start) > slack:
        return False
    if abs(t_end - r_end) > slack:
        return False
    return True


def compute_intron_support(
    read, transcript_intervals: dict[str, Intervals]
) -> list[set[int]]:
    if read["transcript_id"] not in transcript_intervals:
        return list()
    t_intervals = transcript_intervals[read["transcript_id"]]
    support: list[set[int]] = [set() for _ in t_intervals.intervals[:-1]]
    if t_intervals.contig != read["contig"]:
        return support
    for support_set, t_intron in zip(support, generate_introns(t_intervals.intervals)):
        for r_intron in generate_introns(read["intervals"]):
            if close_enough(t_intron, r_intron):
                support_set.add(read.name)
                break
    return support


def aggregate_support(supports):
    L = {len(s) for s in supports}
    assert len(L) == 1
    if L == {0}:
        return 0
    total_support = [set() for _ in supports.iloc[0]]
    for support in supports:
        for S, s in zip(total_support, support):
            S.update(s)
    if min(map(len, total_support)) == 0:
        return 0
    return len(set.union(*total_support))


class CIGAR_OPS_SIMPLE(enum.IntEnum):
    both = 0
    target = 1
    query = 2


op_simply: dict[int, CIGAR_OPS_SIMPLE] = {
    pysam.CSOFT_CLIP: CIGAR_OPS_SIMPLE.query,
    pysam.CINS: CIGAR_OPS_SIMPLE.query,
    pysam.CDEL: CIGAR_OPS_SIMPLE.target,
    pysam.CREF_SKIP: CIGAR_OPS_SIMPLE.target,
    pysam.CMATCH: CIGAR_OPS_SIMPLE.both,
    pysam.CDIFF: CIGAR_OPS_SIMPLE.both,
    pysam.CEQUAL: CIGAR_OPS_SIMPLE.both,
}


@total_ordering
@dataclass
class Interval:
    start: int = 0
    end: int = 0

    def __post_init__(self):
        assert 0 <= self.start <= self.end

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

    def __le__(self, other):
        return (self.start, self.end) <= (other.start, other.end)

    def __gt__(self, other):
        return (self.start, self.end) > (other.start, other.end)

    def __ge__(self, other):
        return (self.start, self.end) >= (other.start, other.end)

    def __len__(self):
        return self.end - self.start


@dataclass
class PairedInterval:
    query: Interval = field(default_factory=Interval)
    target: Interval = field(default_factory=Interval)


@dataclass
class PairedIntervalCigar(PairedInterval):
    cigar: list[tuple[CIGAR_OPS_SIMPLE, int]] = field(default_factory=list)


def canonize_cigar(
    cigartuples: list[tuple[int, int]]
) -> list[tuple[CIGAR_OPS_SIMPLE, int]]:
    simple_cigartuples = [(op_simply[op], l) for op, l in cigartuples]
    canonized_cigar: list[tuple[CIGAR_OPS_SIMPLE, int]] = list()
    for _, g in groupby(
        simple_cigartuples, key=lambda x: x[0] == CIGAR_OPS_SIMPLE.both
    ):
        C: Counter[CIGAR_OPS_SIMPLE] = Counter()
        for op, l in g:
            C[op] += l
        for op, l in sorted(C.items()):
            if l > 0:
                canonized_cigar.append((op, l))
    return canonized_cigar


def get_intervals(
    reference_start: int,
    cigartuples: list[tuple[int, int]],
    cigar_max_del: int = 20,
) -> list[PairedInterval]:
    cigar = canonize_cigar(cigartuples)
    qstart = 0
    qlen = 0
    for op, l in cigar:
        if op in [CIGAR_OPS_SIMPLE.query, CIGAR_OPS_SIMPLE.both]:
            qlen += l

    # list of exonic intervals of the alignment
    intervals: list[PairedIntervalCigar] = list()

    qstart: int = 0  # current interval's start on query
    qend: int = 0  # current interval's end on query
    tstart: int = reference_start  # reference_start is 0-indexed
    tend: int = tstart  # current interval's end on target
    for is_splice, g in groupby(
        cigar,
        key=lambda x: x[0] == CIGAR_OPS_SIMPLE.target and x[1] > cigar_max_del,
    ):
        cur_cigar = list(g)
        for op, l in cur_cigar:
            if op == CIGAR_OPS_SIMPLE.query:
                qend += l
            elif op == CIGAR_OPS_SIMPLE.target:
                tend += l
            elif op == CIGAR_OPS_SIMPLE.both:
                qend += l
                tend += l
        if not is_splice:
            intervals.append(
                PairedIntervalCigar(
                    query=Interval(qstart, qend),
                    target=Interval(tstart, tend),
                    cigar=cur_cigar,
                )
            )
        qstart = qend
        tstart = tend
    final_intervals: list[PairedInterval] = list()
    for interval in intervals:
        qs = interval.query.start
        qe = interval.query.end
        ts = interval.target.start
        te = interval.target.end
        cigar = interval.cigar
        assert qe - qs == (
            S := sum(
                l
                for op, l in cigar
                if op in [CIGAR_OPS_SIMPLE.query, CIGAR_OPS_SIMPLE.both]
            )
        ), (qe - qs, S)
        assert te - ts == (
            S := sum(
                l
                for op, l in cigar
                if op in [CIGAR_OPS_SIMPLE.target, CIGAR_OPS_SIMPLE.both]
            )
        ), (qe - qs, S)
        for op, l in cigar:
            if op == CIGAR_OPS_SIMPLE.both:
                break
            if op == CIGAR_OPS_SIMPLE.query:
                qs += l
            elif op == CIGAR_OPS_SIMPLE.target:
                ts += l
        for op, l in cigar[::-1]:
            if op == CIGAR_OPS_SIMPLE.both:
                break
            if op == CIGAR_OPS_SIMPLE.query:
                qe -= l
            elif op == CIGAR_OPS_SIMPLE.target:
                te -= l
        final_intervals.append(
            PairedInterval(query=Interval(qs, qe), target=Interval(ts, te))
        )
    return final_intervals


def aln_gen(
    bam_path,
) -> Generator[
    tuple[str, str, int, list[tuple[int, int]]],
    None,
    None,
]:
    for aln in tqdm(pysam.AlignmentFile(bam_path).fetch()):
        if aln.is_supplementary or aln.is_secondary or aln.is_unmapped:
            continue
        assert aln.cigartuples is not None
        assert aln.query_name is not None
        assert aln.reference_name is not None

        yield (aln.query_name, aln.reference_name, aln.reference_start, aln.cigartuples)


def process_aln(aln) -> tuple[str, Intervals]:
    query_name, reference_name, reference_start, cigartuples = aln
    return (
        query_name,
        Intervals(
            reference_name,
            [
                (x.target.start, x.target.end)
                for x in get_intervals(reference_start, cigartuples)
            ],
        ),
    )


def get_read_alns(bam_path):
    read_alns: dict[str, Intervals] = dict()
    with multiprocessing.Pool(16) as pool:
        for read_id, aln_intervals in pool.map(process_aln, aln_gen(bam_path)):
            assert read_id not in read_alns
            read_alns[read_id] = aln_intervals
    return read_alns


def get_reads(bam_path, transcript_intervals, truth_tsv):
    read_alns = get_read_alns(bam_path)
    reads = (
        pd.read_csv(
            truth_tsv,
            sep="\t",
        )
        .drop(
            ["cell_barcode"],
            axis=1,
        )
        .fillna("N")
    )
    reads = reads.set_index("read_id", verify_integrity=True)
    reads["contig"] = reads.apply(
        lambda x: read_alns[x.name].contig if x.name in read_alns else "",
        axis=1,
    )
    reads["intervals"] = reads.apply(
        lambda x: read_alns[x.name].intervals if x.name in read_alns else list(),
        axis=1,
    )
    reads["cell_lines"] = reads["cell_lines"].str.split(",")
    reads = reads.explode("cell_lines", ignore_index=True)
    reads["support"] = reads.apply(
        partial(compute_intron_support, transcript_intervals=transcript_intervals),
        axis=1,
    )

    return reads


def main():
    args = parse_args()
    transcript_intervals = get_gtf_data(gtf=args.gtf)
    reads = get_reads(
        bam_path=args.bam,
        transcript_intervals=transcript_intervals,
        truth_tsv=args.truth_tsv,
    )
    transcripts = reads.groupby(["transcript_id", "cell_lines"]).aggregate(
        {"support": aggregate_support}
    )
    transcripts = transcripts[transcripts["support"] >= 3].reset_index()
    transcripts = (
        transcripts.pivot(
            index="transcript_id",
            columns="cell_lines",
            values="support",
        )
        .fillna(0)
        .astype(int)
        .reset_index()
        .rename_axis(None)
        .rename_axis(None, axis=1)
    )
    transcripts["contig"] = transcripts["transcript_id"].apply(
        lambda x: transcript_intervals[x].contig
    )
    transcripts["intervals"] = transcripts["transcript_id"].apply(
        lambda x: ",".join(f"{s}-{e}" for s, e in transcript_intervals[x].intervals)
    )
    transcripts.to_csv(args.o, sep="\t", index=False)


if __name__ == "__main__":
    main()
