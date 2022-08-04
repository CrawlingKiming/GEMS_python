import os
import tokenize

import pandas as pd


def profile_report(stats):
    df = pd.DataFrame(columns=['ncalls', 'time', 'percall', 'tottime', 'totpercall', 'callees'])

    for stat in stats:
        index = "{}:{}({})".format(os.path.basename(stat.code.co_filename), stat.code.co_firstlineno, stat.code.co_name)
        df.loc[index, 'ncalls'] = stat.callcount
        df.loc[index, 'time'] = stat.inlinetime
        df.loc[index, 'percall'] = stat.inlinetime / stat.callcount
        df.loc[index, 'tottime'] = stat.totaltime
        df.loc[index, 'totpercall'] = stat.totaltime / stat.callcount
        if stat.calls is not None:
            df.loc[index, 'callees'] = ", ".join([x.code.co_name for x in stat.calls])

    df.index.name = "filename:lineno(function)"

    df.sort_values('time', axis=0, ascending=False, inplace=True)
    df.to_excel("profile_result.xlsx", engine="openpyxl")


def line_profile_report(lstats, output_unit=0.001):
    scalar = lstats.unit / output_unit

    stats = lstats.timings

    writer = pd.ExcelWriter("line_profile_result.xlsx", engine="openpyxl")

    df_summary = pd.DataFrame(columns=['total_time'])

    for (filename, lineno, name), timings in sorted(stats.items()):

        if len(timings) == 0:
            continue

        start_lineno = lineno - 1
        last_lineno = timings[-1][0]

        df = pd.DataFrame(columns=['code', 'hits', 'time'])

        with tokenize.open(filename) as fp:
            lines = fp.readlines()

        for i in range(start_lineno, last_lineno, 1):
            line = lines[i].replace("\n", "")
            df.loc[i + 1, 'code'] = line

        for lineno, nhits, time in timings:
            df.loc[lineno, 'hits'] = nhits
            df.loc[lineno, 'time'] = time * scalar

        df['per_hits'] = df['time'] / df['hits']
        df['ratio'] = df['time'] / df['time'].sum()
        df.index.name = 'lineno'

        if name in df_summary.index:
            name += "0"

        df.to_excel(writer, sheet_name=name)
        df_summary.loc[name, 'total_time'] = df['time'].sum()

    df_summary.index.name = 'function'
    df_summary.to_excel(writer, sheet_name="summary")
    writer.save()
