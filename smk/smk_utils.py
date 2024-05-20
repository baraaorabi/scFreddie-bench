import datetime
import typing
from tqdm import tqdm
import numpy as np

np.random.seed(42)

time_header_str = ",".join(
    [
        "stamp",
        "sample",
        "proc",
        "real_min",
        "user_min",
        "sys_min",
        "mem_GB",
        "threads",
    ]
)


def format_gnu_time_string(
    process="",
    sample="{wildcards.sample}",
    threads="{threads}",
):
    fields = list()
    fields.append(
        (f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", "%s", "")
    )
    fields.append((f"{sample}", "%s", ""))
    fields.append((f"{process}", "%s", ""))
    fields.append(("%e", "%.2f", "/60"))
    fields.append(("%U", "%.2f", "/60"))
    fields.append(("%S", "%.2f", "/60"))
    fields.append(("%M", "%.2f", "/(1024*1024)"))
    fields.append((f"{threads}", "%d", ""))

    time_format = ",".join([x[0] for x in fields])
    printf_format = ",".join([x[1] for x in fields]) + "\\n"
    printf_args = ",".join([f"${i}{x[2]}" for i, x in enumerate(fields, start=1)])

    awk_cmd = (
        'awk \'BEGIN{{FS=","}} {{printf "' + printf_format + '",' + printf_args + "}}'"
    )
    return f'$(which time) -f "{time_format}" -o >({awk_cmd} >> {{input.time}}) '


def sample_gtf(gtf: str, rate: float, outfile: typing.TextIO):
    assert 0.0 <= rate <= 1.0
    genes: dict[str, tuple[str, dict[str, list[str]]]] = dict()
    tid_to_gid: dict[str, str] = dict()
    all_gids = list()
    for line in tqdm(open(gtf)):
        line = line.rstrip("\n")
        if line.startswith("#"):
            print(line, file=outfile)
            continue
        fields = line.split("\t")
        info = fields[8]
        info = [x.strip().split(" ") for x in info.strip(";").split(";")]
        info = {x[0]: x[1].strip('"') for x in info}
        if fields[2] == "gene":
            genes[info["gene_id"]] = (line, dict())
            all_gids.append(info["gene_id"])
            continue
        if fields[2] == "transcript":
            tid_to_gid[info["transcript_id"]] = info["gene_id"]
            genes[info["gene_id"]][1][info["transcript_id"]] = list()
        genes[info["gene_id"]][1][info["transcript_id"]].append(line)

    N = len(tid_to_gid)
    K = round(N * rate)
    print(f"Sampling {K} transcripts from {N} transcripts ({rate:.2%})")
    sample_tids = np.random.choice(
        list(tid_to_gid.keys()),
        size=round(len(tid_to_gid) * rate),
        replace=False,
    )

    sample_gids = {tid_to_gid[tid] for tid in sample_tids}
    for gid in all_gids:
        if not gid in sample_gids:
            continue
        print(genes[gid][0], file=outfile)
        for tid, lines in genes[gid][1].items():
            if not tid in sample_tids:
                continue
            for line in lines:
                print(line, file=outfile)
