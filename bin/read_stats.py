#!/usr/bin/env python3
"""
Read-level QC summary (read length + mean Q-score) for one or more FASTQ
stages -- e.g. before vs. after chopper filtering -- in the spirit of NanoPlot.

Pure standard library (no pandas/numpy/matplotlib), so it runs in any python3
container and stand-alone on a laptop. Reads .fastq / .fastq.gz (detected by
magic bytes, like `zcat -f`).

Two subcommands:

  summarize  FASTQ(s) -> per-stage stats TSV + histogram TSV (+ optional HTML)
  report     stats/hist TSVs from many samples -> one combined HTML report

Stand-alone, one sample, before/after, one command:

  read_stats.py summarize --sample bc01 \\
      --stage raw=bc01.fastq.gz --stage filtered=bc01.filtered.fastq.gz \\
      --out-prefix bc01 --html bc01.read_qc.html

Many samples (what the pipeline does): summarize each, then combine:

  read_stats.py report --stats *.read_stats.tsv --hist *.read_hist.tsv \\
      --out-table read_stats.tsv --out-hist read_hist.tsv --out-html read_qc.html

Stage names are free-form (raw, filtered, trimmed, ...); order on the command
line is the order shown in tables and plots, and the first stage is the
baseline for "% retained".

Per-read Q is the mean *error probability* converted back to Phred
(-10*log10(mean(10^(-q/10)))), the same convention chopper and NanoPlot use --
not the plain average of the per-base Q values.
"""
import argparse
import csv
import gzip
import html
import math
import sys
from array import array
from collections import Counter, defaultdict

STATS_COLS = [
    "sample", "stage", "n_reads", "total_bases", "mean_len", "median_len",
    "n50", "min_len", "max_len", "mean_q", "median_q",
    "pct_q10", "pct_q15", "pct_q20",
]
HIST_COLS = ["sample", "stage", "metric", "bin_start", "bin_end", "reads"]
POOLED = "All samples"

# per-byte error probability, indexed by the raw ASCII value of a quality char
_ERR = [10 ** (-max(b - 33, 0) / 10) for b in range(256)]


# --------------------------------------------------------------------------
# FASTQ parsing + per-stage summary
# --------------------------------------------------------------------------
def _open_maybe_gz(path):
    fh = open(path, "rb")
    magic = fh.read(2)
    fh.seek(0)
    if magic == b"\x1f\x8b":
        return gzip.open(fh, "rb")
    return fh


def read_fastq(path):
    """Yield (length, quality_bytes) per record; assumes 4-line FASTQ, which
    is what ONT basecallers write."""
    with _open_maybe_gz(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            if not header.strip():
                continue
            if not header.startswith(b"@"):
                raise ValueError(f"{path}: expected '@' header line, got {header[:30]!r}")
            seq = fh.readline()
            fh.readline()
            qual = fh.readline()
            if not qual:
                raise ValueError(f"{path}: truncated record after {header.strip()!r}")
            yield len(seq.rstrip(b"\r\n")), qual.rstrip(b"\r\n")


def mean_q(qual):
    return -10 * math.log10(max(sum(map(_ERR.__getitem__, qual)) / len(qual), 1e-10))


def median(sorted_vals):
    n = len(sorted_vals)
    if n == 0:
        return None
    mid = n // 2
    return sorted_vals[mid] if n % 2 else (sorted_vals[mid - 1] + sorted_vals[mid]) / 2


def n50(lengths):
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= total:
            return length
    return None


def summarize_fastq(path):
    lengths = array("I")
    quals = array("d")
    for length, qual in read_fastq(path):
        lengths.append(length)
        if length:
            quals.append(mean_q(qual))
    return lengths, quals


def stats_row(sample, stage, lengths, quals):
    n = len(lengths)
    sl, sq = sorted(lengths), sorted(quals)
    frac = lambda thr: (100 * sum(q >= thr for q in quals) / len(quals)) if quals else None
    return {
        "sample": sample, "stage": stage, "n_reads": n,
        "total_bases": sum(lengths),
        "mean_len": sum(lengths) / n if n else None,
        "median_len": median(sl), "n50": n50(lengths) if n else None,
        "min_len": sl[0] if n else None, "max_len": sl[-1] if n else None,
        "mean_q": sum(quals) / len(quals) if quals else None,
        "median_q": median(sq),
        "pct_q10": frac(10), "pct_q15": frac(15), "pct_q20": frac(20),
    }


def hist_rows(sample, stage, lengths, quals, len_bin, q_bin):
    rows = []
    for metric, vals, width in (("length", lengths, len_bin), ("qscore", quals, q_bin)):
        counts = Counter(int(v // width) for v in vals)
        for idx in sorted(counts):
            rows.append({
                "sample": sample, "stage": stage, "metric": metric,
                "bin_start": round(idx * width, 6),
                "bin_end": round((idx + 1) * width, 6),
                "reads": counts[idx],
            })
    return rows


# --------------------------------------------------------------------------
# TSV io
# --------------------------------------------------------------------------
def _fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.2f}"
    return str(v)


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([_fmt(r[c]) for c in cols])


def _num(v):
    if v == "NA":
        return None
    try:
        return int(v)
    except ValueError:
        return float(v)


def read_tsv(paths, numeric):
    rows = []
    for p in paths:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for c in numeric:
                    r[c] = _num(r[c])
                rows.append(r)
    return rows


STATS_NUMERIC = STATS_COLS[2:]
HIST_NUMERIC = ["bin_start", "bin_end", "reads"]


# --------------------------------------------------------------------------
# HTML report (inline SVG -- no JS, no CDN, works offline / as a GCS file)
# --------------------------------------------------------------------------
def _si(n):
    if n is None:
        return "NA"
    for div, suffix in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(n) >= div:
            return f"{n / div:.3g}{suffix}"
    return f"{n:.3g}"


def _ticks(hi, target=5):
    if hi <= 0:
        return [0]
    raw = hi / target
    mag = 10 ** math.floor(math.log10(raw))
    step = next(m * mag for m in (1, 2, 5, 10) if m * mag >= raw)
    return [i * step for i in range(int(hi / step) + 1)]


def svg_hist(series, x_label, x_max=None, w=460, h=230):
    """series: list of (stage_index, [(bin_start, bin_end, reads), ...]).
    Overlaid step outlines, one per stage, so before/after is directly
    comparable on a shared axis."""
    pad_l, pad_r, pad_t, pad_b = 52, 12, 10, 38
    bins = [b for _, rows in series for b in rows]
    if not bins:
        return '<p class="muted">No reads.</p>'
    x_lo = min(b[0] for b in bins)
    x_hi = x_max if x_max is not None else max(b[1] for b in bins)
    if x_hi <= x_lo:
        x_hi = x_lo + 1
    y_hi = max(b[2] for b in bins if b[0] < x_hi)
    y_ticks = _ticks(y_hi)
    y_top = y_ticks[-1] if y_ticks[-1] >= y_hi else y_ticks[-1] + (y_ticks[1] - y_ticks[0])
    pw, ph = w - pad_l - pad_r, h - pad_t - pad_b
    X = lambda v: pad_l + (min(v, x_hi) - x_lo) / (x_hi - x_lo) * pw
    Y = lambda v: pad_t + ph - v / y_top * ph

    out = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{html.escape(x_label)} histogram">']
    for t in _ticks(y_top):
        out.append(f'<line class="grid" x1="{pad_l}" x2="{w - pad_r}" y1="{Y(t):.1f}" y2="{Y(t):.1f}"/>'
                   f'<text class="tick" x="{pad_l - 6}" y="{Y(t) + 3:.1f}" text-anchor="end">{_si(t)}</text>')
    for t in _ticks(x_hi):
        if t >= x_lo:
            out.append(f'<text class="tick" x="{X(t):.1f}" y="{h - pad_b + 14}" text-anchor="middle">{_si(t)}</text>')
    for idx, rows in series:
        pts = []
        for start, end, reads in sorted(rows):
            if start >= x_hi:
                break
            pts += [f"{X(start):.1f},{Y(reads):.1f}", f"{X(end):.1f},{Y(reads):.1f}"]
        if not pts:
            continue
        base = f"{pad_t + ph:.1f}"
        first_x, last_x = pts[0].split(",")[0], pts[-1].split(",")[0]
        out.append(f'<polygon class="area s{idx % 4}" points="{first_x},{base} {" ".join(pts)} {last_x},{base}"/>'
                   f'<polyline class="line s{idx % 4}" points="{" ".join(pts)}"/>')
    out.append(f'<line class="axis" x1="{pad_l}" x2="{w - pad_r}" y1="{pad_t + ph}" y2="{pad_t + ph}"/>'
               f'<line class="axis" x1="{pad_l}" x2="{pad_l}" y1="{pad_t}" y2="{pad_t + ph}"/>'
               f'<text class="label" x="{pad_l + pw / 2}" y="{h - 4}" text-anchor="middle">{html.escape(x_label)}</text>'
               f'<text class="label" transform="translate(12 {pad_t + ph / 2}) rotate(-90)" text-anchor="middle">reads</text>')
    out.append("</svg>")
    return "".join(out)


def _x_max_len(series):
    """Clip the length axis at the 99th percentile of the widest stage so a
    few very long reads don't squash the amplicon-sized bulk into one bar."""
    hi = 0
    for _, rows in series:
        total = sum(r[2] for r in rows)
        running = 0
        for start, end, reads in sorted(rows):
            running += reads
            if running >= 0.99 * total:
                hi = max(hi, end)
                break
    return hi or None


CSS = """
:root{--bg:#fff;--fg:#1c2024;--muted:#667;--border:#dde1e6;--card:#f7f8fa;
--s0:#2b6cb0;--s1:#d9822b;--s2:#2f9e6e;--s3:#8a5bc7}
@media(prefers-color-scheme:dark){:root{--bg:#15181c;--fg:#e6e8ea;--muted:#9aa4ae;
--border:#2c333a;--card:#1c2025;--s0:#63a4e8;--s1:#f0a45a;--s2:#5cc79a;--s3:#b593e6}}
body{background:var(--bg);color:var(--fg);font:14px/1.5 system-ui,sans-serif;margin:0 auto;max-width:1000px;padding:16px}
h1{font-size:20px}h2{font-size:16px;margin:0}
.muted{color:var(--muted)}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
th,td{padding:4px 8px;border-bottom:1px solid var(--border);text-align:right;white-space:nowrap}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}
.scroll{overflow-x:auto}
details{border:1px solid var(--border);border-radius:6px;background:var(--card);margin:12px 0;padding:8px 12px}
summary{cursor:pointer}
.plots{display:flex;flex-wrap:wrap;gap:12px;margin-top:8px}
.plots svg{flex:1 1 320px;max-width:100%;height:auto}
.legend span{display:inline-block;margin-right:14px}
.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px}
svg .grid{stroke:var(--border)}svg .axis{stroke:var(--muted)}
svg .tick,svg .label{fill:var(--muted);font-size:10px}
svg .area{opacity:.22}svg .line{fill:none;stroke-width:1.5}
.s0{fill:var(--s0);stroke:var(--s0)}.s1{fill:var(--s1);stroke:var(--s1)}
.s2{fill:var(--s2);stroke:var(--s2)}.s3{fill:var(--s3);stroke:var(--s3)}
polyline.line{fill:none!important}
.legend i.s0{background:var(--s0)}.legend i.s1{background:var(--s1)}
.legend i.s2{background:var(--s2)}.legend i.s3{background:var(--s3)}
"""


def _stats_table(rows, stages):
    base = {}
    for r in rows:
        base.setdefault(r["sample"], r)
    head = ["Sample", "Stage", "Reads", "% reads kept", "Bases", "% bases kept", "Mean len",
            "Median len", "N50", "Mean Q", "Median Q", "% Q≥10", "% Q≥15", "% Q≥20"]
    out = ['<div class="scroll"><table><thead><tr>'
           + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead><tbody>"]
    for r in rows:
        b = base[r["sample"]]
        pct = lambda a, z: "NA" if not z else f"{100 * a / z:.1f}"
        cells = [
            html.escape(r["sample"]), html.escape(r["stage"]), f'{r["n_reads"]:,}',
            pct(r["n_reads"], b["n_reads"]), f'{r["total_bases"]:,}',
            pct(r["total_bases"], b["total_bases"]),
            _fmt(r["mean_len"]), _fmt(r["median_len"]), _fmt(r["n50"]),
            _fmt(r["mean_q"]), _fmt(r["median_q"]),
            _fmt(r["pct_q10"]), _fmt(r["pct_q15"]), _fmt(r["pct_q20"]),
        ]
        out.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def build_html(stats, hists, title="Read QC summary"):
    stages = []
    for r in stats:
        if r["stage"] not in stages:
            stages.append(r["stage"])
    samples = sorted({r["sample"] for r in stats})

    # sample -> stage -> metric -> [(start, end, reads)]
    by = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    pooled = defaultdict(lambda: defaultdict(Counter))
    for r in hists:
        by[r["sample"]][r["stage"]][r["metric"]].append((r["bin_start"], r["bin_end"], r["reads"]))
        pooled[r["stage"]][r["metric"]][(r["bin_start"], r["bin_end"])] += r["reads"]
    for stage, metrics in pooled.items():
        for metric, counts in metrics.items():
            by[POOLED][stage][metric] = [(s, e, n) for (s, e), n in counts.items()]

    legend = '<div class="legend">' + "".join(
        f'<span><i class="s{i % 4}"></i>{html.escape(s)}</span>' for i, s in enumerate(stages)) + "</div>"

    def section(name, open_):
        series = lambda metric: [(i, by[name][s][metric]) for i, s in enumerate(stages) if by[name][s][metric]]
        len_s = series("length")
        return (f'<details{" open" if open_ else ""}><summary><h2 style="display:inline">{html.escape(name)}</h2></summary>'
                f'{legend}<div class="plots">'
                f'{svg_hist(len_s, "read length (bp)", _x_max_len(len_s))}'
                f'{svg_hist(series("qscore"), "mean read Q-score")}</div></details>')

    ordered = sorted(stats, key=lambda r: (r["sample"], stages.index(r["stage"])))
    pooled_rows = []
    for stage in stages:
        rs = [r for r in stats if r["stage"] == stage]
        pooled_rows.append({c: None for c in STATS_COLS} | {
            "sample": POOLED, "stage": stage,
            "n_reads": sum(r["n_reads"] for r in rs),
            "total_bases": sum(r["total_bases"] for r in rs),
        })
        pooled_rows[-1]["mean_len"] = (pooled_rows[-1]["total_bases"] / pooled_rows[-1]["n_reads"]
                                       if pooled_rows[-1]["n_reads"] else None)
    body = [f"<h1>{html.escape(title)}</h1>",
            '<p class="muted">Per-read length and mean Q-score at each stage; the first stage is the '
            "baseline for “% kept”. Q is the mean error probability converted back to Phred.</p>",
            _stats_table(pooled_rows + ordered, stages)]
    if len(samples) > 1:
        body.append(section(POOLED, True))
    body += [section(s, len(samples) == 1) for s in samples]
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>{html.escape(title)}</title><style>{CSS}</style></head><body>"
            + "".join(body) + "</body></html>")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def cmd_summarize(args):
    stats, hists = [], []
    for spec in args.stage:
        stage, sep, path = spec.partition("=")
        if not sep or not stage or not path:
            sys.exit(f"--stage expects NAME=FASTQ, got {spec!r}")
        lengths, quals = summarize_fastq(path)
        stats.append(stats_row(args.sample, stage, lengths, quals))
        hists += hist_rows(args.sample, stage, lengths, quals, args.len_bin, args.q_bin)
    write_tsv(f"{args.out_prefix}.read_stats.tsv", STATS_COLS, stats)
    write_tsv(f"{args.out_prefix}.read_hist.tsv", HIST_COLS, hists)
    if args.html:
        with open(args.html, "w") as fh:
            fh.write(build_html(stats, hists, args.title))


def cmd_report(args):
    stats = read_tsv(sorted(args.stats), STATS_NUMERIC)
    hists = read_tsv(sorted(args.hist), HIST_NUMERIC)
    if not stats:
        sys.exit("no stats rows found in --stats inputs")
    write_tsv(args.out_table, STATS_COLS, stats)
    write_tsv(args.out_hist, HIST_COLS, hists)
    with open(args.out_html, "w") as fh:
        fh.write(build_html(stats, hists, args.title))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("summarize", help="FASTQ stage(s) -> stats + histogram TSVs")
    s.add_argument("--sample", required=True)
    s.add_argument("--stage", action="append", required=True, metavar="NAME=FASTQ",
                   help="repeatable, in display order; first is the baseline")
    s.add_argument("--out-prefix", required=True, help="writes PREFIX.read_stats.tsv and PREFIX.read_hist.tsv")
    s.add_argument("--len-bin", type=float, default=10, help="length histogram bin width, bp (default 10)")
    s.add_argument("--q-bin", type=float, default=0.5, help="Q-score histogram bin width (default 0.5)")
    s.add_argument("--html", help="also write a stand-alone HTML report for this sample")
    s.add_argument("--title", default="Read QC summary")
    s.set_defaults(func=cmd_summarize)

    r = sub.add_parser("report", help="combine per-sample TSVs into one report")
    r.add_argument("--stats", nargs="+", required=True)
    r.add_argument("--hist", nargs="+", required=True)
    r.add_argument("--out-table", default="read_stats.tsv")
    r.add_argument("--out-hist", default="read_hist.tsv")
    r.add_argument("--out-html", default="read_qc_summary.html")
    r.add_argument("--title", default="Read QC summary")
    r.set_defaults(func=cmd_report)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
