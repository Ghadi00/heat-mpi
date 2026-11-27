#!/usr/bin/env python3
"""fit_pingpong.py -- fit t_msg(m) = t_s + m * t_w to ping-pong results.

Usage:
    python3 scripts/fit_pingpong.py results/pingpong/<run-folder>

For every *.csv file in the folder (one per pair of hosts) this prints
  * t_s  : latency, the mean one-way time of the small messages (m <= 1 KiB)
  * t_w  : time per byte, the least-squares slope of the one-way time over
           the large messages (m >= 1 MiB); bandwidth = 1 / t_w
The one-way time is half the MEDIAN round trip by default, so that a single
delayed round trip (a Wi-Fi hiccup, a relay switch) does not distort the
fit; --stat mean uses half the average round trip instead.
and writes
  * <folder>/summary.md   a table for the report
  * <folder>/pingpong.png one-way time vs message size, measured and model

Following Section 8.4 of the background notes: the latency comes from small
messages and the bandwidth from large ones.
"""
import glob
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SMALL_MAX = 1024          # bytes: messages used for t_s
LARGE_MIN = 1024 * 1024   # bytes: messages used for t_w


def read_run(path):
    with open(path) as f:
        first = f.readline().strip()
    label = first.replace("# pingpong", "").strip() if first.startswith("#") else os.path.basename(path)
    with open(path) as f:
        rows = [ln for ln in f if ln.strip() and not ln.startswith("#")]
    data = np.genfromtxt(rows, delimiter=",", names=True)
    return label, np.atleast_1d(data)


def fit(data, stat="median"):
    m = data["bytes"].astype(float)
    t = (data["rtt_median_s"].astype(float) / 2.0) if stat == "median" \
        else data["oneway_avg_s"].astype(float)
    small = m <= SMALL_MAX
    large = m >= LARGE_MIN
    ts = float(np.mean(t[small])) if small.any() else float("nan")
    if large.sum() >= 2:
        tw, _ = np.polyfit(m[large], t[large], 1)
    else:
        tw = float("nan")
    return ts, float(tw)


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--stat")]
    stat = "mean" if "--stat=mean" in sys.argv or ("--stat" in sys.argv and "mean" in sys.argv) else "median"
    args = [x for x in args if x not in ("mean", "median")]
    if len(args) != 1:
        print(__doc__)
        sys.exit(1)
    folder = args[0]
    files = sorted(glob.glob(os.path.join(folder, "*.csv")))
    if not files:
        sys.exit(f"no .csv files in {folder}")

    rows = []
    fig, ax = plt.subplots(figsize=(7, 5))
    for path in files:
        label, data = read_run(path)
        if data.size == 0:
            print(f"skipping empty file {path}")
            continue
        ts, tw = fit(data, stat)
        bw = 1.0 / tw / 1e6 if tw > 0 else float("nan")
        rows.append((label, ts, tw, bw))

        m = data["bytes"].astype(float)
        y = data["rtt_median_s"] / 2.0 if stat == "median" else data["oneway_avg_s"]
        line, = ax.loglog(m, y, "o", ms=4, label=label)
        ax.loglog(m, ts + m * tw, "-", color=line.get_color(), lw=1, alpha=0.7)

    ax.set_xlabel("message size m [bytes]")
    ax.set_ylabel("one-way time [s]")
    ax.set_title("Ping-pong: measured (dots) and model t_s + m t_w (lines)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    png = os.path.join(folder, "pingpong.png")
    fig.savefig(png, dpi=150)

    lines = [
        "| Pair | t_s (latency) | t_w (s/byte) | Bandwidth 1/t_w |",
        "|---|---|---|---|",
    ]
    for label, ts, tw, bw in rows:
        lines.append(f"| {label} | {ts * 1e6:.1f} µs | {tw:.3e} | {bw:.2f} MB/s |")
    table = "\n".join(lines)
    print(table)
    with open(os.path.join(folder, "summary.md"), "w") as f:
        f.write(table + "\n\nAdd a column by hand: direct or DERP (see the *.path.txt files).\n")
    print(f"\nwrote {png} and {os.path.join(folder, 'summary.md')}")


if __name__ == "__main__":
    main()
