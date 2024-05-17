import datetime


header_str = ",".join(
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
