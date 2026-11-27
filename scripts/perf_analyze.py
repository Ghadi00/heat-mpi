#!/usr/bin/env python3
"""perf_analyze.py -- tables, plots and the performance model for Part 4.

Usage:
    python3 scripts/perf_analyze.py results/perf/<folder> \
        --pingpong results/pingpong/<folder>        # t_s, t_w measured in Part 0
    # or give the network parameters directly:
    python3 scripts/perf_analyze.py results/perf/<folder> --ts 0.07 --tw 4e-6

Reads <folder>/runs.csv (written by perf_run.py), takes the MEDIAN of the repetitions
of every configuration (and the min/max as the spread) and writes into <folder>/analysis/:

  summary.md                all tables: compute rate, strong and weak scaling with spread,
                            communication share, model parameters, model vs measurement,
                            and the model's prediction for a real cluster
  baseline.png              time per node per step on each VM (compute rate)
  strong_N<N>.png           time per step, speed-up S = T1/TP, efficiency E = S/P,
                            measured (markers, min-max bars) and model (dashed)
  weak.png                  weak scaling: time per step vs P, measured and model
  comm_share.png            fraction of the time spent in the halo exchange
  real_cluster.png          model speed-up with t_s = 1 us, 10 GB/s (Section 8.3 question)

THE MODEL (background notes, Task 8.2; time per step, slowest rank decides):
  compute   t_node * (number of interior nodes of the biggest piece)
            t_node = measured serial time per node per step, from the baseline
            (the SLOWEST VM that takes part, since every rank waits for the slowest)
  strips    + 2 * (t_s + 8 (N+1) t_w)                  for P > 1
            two Sendrecv calls per step, each moves one row of N+1 doubles
            (send and receive overlap, so one message time each)
  blocks    + 2 * (t_s + 8 nx t_w)   if Py > 1   (rows, nx = widest block)
            + 2 * (t_s + 8 ny t_w)   if Px > 1   (columns, ny = tallest block)
  t_s, t_w: the cross-VM values from the ping-pong fit as soon as any neighbouring
            ranks are on different VMs (the run used more than one VM: always with
            --map-by node, with --map-by slot once P exceeds the first VM's slots),
            the same-VM values otherwise.
            With several VM pairs the WORST link is used (--link mean to average).
"""
import argparse
import csv
import glob
import math
import os
import statistics
import sys
from collections import defaultdict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ------------------------------------------------------------------ data --
def load(folder):
    rows = []
    with open(os.path.join(folder, "runs.csv")) as f:
        for r in csv.DictReader(f):
            if r["ok"] != "1":
                continue
            for k in ("N", "P", "steps", "rep"):
                r[k] = int(r[k])
            for k in ("t_total", "t_comp", "t_comm"):
                r[k] = float(r[k])
            rows.append(r)
    return rows


def group(rows, keys):
    g = defaultdict(list)
    for r in rows:
        g[tuple(r[k] for k in keys)].append(r)
    return g


def stats(runs, field="t_total", per_step=True):
    xs = [r[field] / (r["steps"] if per_step else 1) for r in runs]
    return statistics.median(xs), min(xs), max(xs), len(xs)


def split_count(k, N, parts):
    return (N - 1) // parts + (1 if k < (N - 1) % parts else 0)


def dims_create(P):
    """Same result as MPI_Dims_create(P, 2): the most square factorisation, dims[0] >= dims[1]."""
    best = (P, 1)
    for py in range(1, int(math.isqrt(P)) + 1):
        if P % py == 0:
            best = (P // py, py)
    return best


# ----------------------------------------------------------------- model --
class Model:
    def __init__(self, tnode_by_host_N, net, local, slots_first, link="worst"):
        self.tnode = tnode_by_host_N     # {(host, N): seconds per node per step}
        self.net = net                   # (ts, tw) across VMs
        self.local = local               # (ts, tw) within one VM
        self.slots_first = slots_first

    def t_node(self, N, hosts):
        vals = [v for (h, n), v in self.tnode.items() if n == N and (not hosts or h in hosts)]
        if not vals:   # nearest N measured
            Ns = sorted({n for (_, n) in self.tnode})
            if not Ns:
                return float("nan")
            n2 = min(Ns, key=lambda n: abs(math.log(n / N)))
            vals = [v for (h, n), v in self.tnode.items() if n == n2 and (not hosts or h in hosts)] or \
                   [v for (h, n), v in self.tnode.items() if n == n2]
        return max(vals)

    def step(self, decomp, mapping, N, P, hosts=None, ts_tw=None):
        tn = self.t_node(N, hosts)
        if decomp == "serial" or P == 1:
            return tn * (N - 1) ** 2, 0.0
        if ts_tw is None:
            # across VMs as soon as the run uses more than one VM; if that is
            # unknown, decide from the placement rule
            cross = (len(hosts) > 1) if hosts else (mapping == "node" or P > self.slots_first)
            ts_tw = self.net if cross else self.local
        ts, tw = ts_tw
        if decomp == "strips":
            comp = tn * (N - 1) * split_count(0, N, P)
            comm = 2 * (ts + 8 * (N + 1) * tw)
        else:
            Px, Py = dims_create(P)
            nx, ny = split_count(0, N, Px), split_count(0, N, Py)
            comp = tn * nx * ny
            comm = (2 * (ts + 8 * nx * tw) if Py > 1 else 0) + (2 * (ts + 8 * ny * tw) if Px > 1 else 0)
        return comp, comm


def pingpong_params(folder, link):
    from fit_pingpong import read_run, fit
    cross, same = [], []
    for path in sorted(glob.glob(os.path.join(folder, "*.csv"))):
        label, data = read_run(path)
        if data.size == 0:
            continue
        ts, tw = fit(data)
        a, _, b = label.partition("<->")
        (same if a.strip() == b.strip() else cross).append((ts, tw, label))
    if not cross:
        sys.exit(f"no cross-VM ping-pong files in {folder}")
    pick = (lambda xs: (max(x[0] for x in xs), max(x[1] for x in xs))) if link == "worst" else \
           (lambda xs: (statistics.mean(x[0] for x in xs), statistics.mean(x[1] for x in xs)))
    net = pick(cross)
    loc = (min(x[0] for x in same), min(x[1] for x in same)) if same else (5e-7, 1.5e-10)
    return net, loc, cross, same


# ---------------------------------------------------------------- output --
def fmt_t(x):
    if x != x:
        return "—"
    if x >= 1:
        return f"{x:.3f} s"
    if x >= 1e-3:
        return f"{x * 1e3:.2f} ms"
    return f"{x * 1e6:.1f} µs"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--pingpong", help="ping-pong results folder (Part 0) for t_s and t_w")
    ap.add_argument("--ts", type=float, help="cross-VM latency [s] (instead of --pingpong)")
    ap.add_argument("--tw", type=float, help="cross-VM time per byte [s/B]")
    ap.add_argument("--link", choices=["worst", "mean"], default="worst")
    ap.add_argument("--hosts", default=None, help="hostfile used (default: ./hosts)")
    ap.add_argument("--real-ts", type=float, default=1e-6)
    ap.add_argument("--real-bw", type=float, default=10e9, help="bytes/s for the real-cluster prediction")
    a = ap.parse_args()

    rows = load(a.folder)
    if not rows:
        sys.exit("no successful runs in runs.csv")
    out = os.path.join(a.folder, "analysis")
    os.makedirs(out, exist_ok=True)
    md = []

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    hostfile = a.hosts or os.path.join(root, "hosts")
    slots_first = 2
    if os.path.exists(hostfile):
        from perf_run import read_hosts
        slots_first = read_hosts(hostfile)[0][1]

    # ---- baseline: compute rate per VM
    base = group([r for r in rows if r["exp"] == "baseline"], ["host", "N"])
    tnode = {}
    for (h, N), rs in base.items():
        med = statistics.median(r["t_total"] / (r["steps"] * (N - 1) ** 2) for r in rs)
        tnode[(h, N)] = med
    # fall back on the strong-scaling P=1 runs if there is no baseline
    for (h, N), rs in group([r for r in rows if r["decomp"] == "serial" and r["exp"] != "baseline"],
                            ["host", "N"]).items():
        tnode.setdefault((h, N), statistics.median(r["t_total"] / (r["steps"] * (N - 1) ** 2) for r in rs))

    if base:
        hosts = sorted({h for h, _ in base})
        Ns = sorted({N for _, N in base})
        md += ["## 1. Single-VM baseline (compute rate)", "",
               "Serial time per step and per node-update; median of the repetitions, spread in brackets.", "",
               "| N | " + " | ".join(f"{h}: per step | {h}: ns/node" for h in hosts) + " |",
               "|---|" + "---|---|" * len(hosts)]
        for N in Ns:
            cells = []
            for h in hosts:
                rs = base.get((h, N))
                if not rs:
                    cells += ["—", "—"]
                    continue
                m, lo, hi, n = stats(rs)
                cells += [f"{fmt_t(m)} ({fmt_t(lo)}–{fmt_t(hi)})", f"{m / (N - 1) ** 2 * 1e9:.2f}"]
            md.append(f"| {N} | " + " | ".join(cells) + " |")
        md.append("")
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for h in hosts:
            xs = [N for N in Ns if (h, N) in tnode]
            ax.plot(xs, [tnode[(h, N)] * 1e9 for N in xs], "o-", label=h)
        ax.set_xscale("log", base=2)
        ax.set_xlabel("N")
        ax.set_ylabel("ns per node per step")
        ax.set_title("Serial compute rate on each VM")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out, "baseline.png"), dpi=150)
        plt.close(fig)

    # ---- model parameters
    model = None
    if a.pingpong or (a.ts is not None and a.tw is not None):
        if a.pingpong:
            net, loc, cross, same = pingpong_params(a.pingpong, a.link)
            src = f"ping-pong folder `{a.pingpong}`, {a.link} cross-VM pair"
        else:
            net, loc, src = (a.ts, a.tw), (5e-7, 1.5e-10), "command line"
        model = Model(tnode, net, loc, slots_first, a.link)
        md += ["## Model parameters", "",
               f"* t_s, t_w across VMs: {fmt_t(net[0])}, {net[1]:.3e} s/B "
               f"(bandwidth {1 / net[1] / 1e6:.3f} MB/s), from {src}",
               f"* t_s, t_w within a VM: {fmt_t(loc[0])}, {loc[1]:.3e} s/B",
               "* t_node: median serial time per node per step from the baseline, slowest VM taking part", ""]

    # ---- strong scaling
    strong = [r for r in rows if r["exp"] == "strong"]
    if strong:
        md += ["## 2. Strong scaling", "",
               "Time per step (median, spread min–max), speed-up S = T1/TP, efficiency E = S/P, "
               "communication share t_comm/t_total. T1 = serial run on the first host.", ""]
        g = group(strong, ["decomp", "map", "N", "P"])
        for N in sorted({r["N"] for r in strong}):
            t1runs = g.get(("serial", "-", N, 1))
            if not t1runs:
                continue
            T1 = stats(t1runs)[0]
            md += [f"### N = {N}  (T1 = {fmt_t(T1)} per step)", "",
                   "| decomp | map | P | time/step (spread) | S | E | comm share |" +
                   (" model time/step | model S |" if model else ""),
                   "|---|---|---|---|---|---|---|" + ("---|---|" if model else "")]
            fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
            series = sorted({(r["decomp"], r["map"]) for r in strong if r["decomp"] != "serial"})
            colors = {("strips", "node"): "C0", ("strips", "slot"): "C1",
                      ("blocks", "node"): "C2", ("blocks", "slot"): "C3"}
            for d, m in series:
                Ps = sorted(P for (dd, mm, NN, P) in g if dd == d and mm == m and NN == N)
                Ps = [1] + Ps
                med, lo, hi, S, E = [], [], [], [], []
                for P in Ps:
                    rs = t1runs if P == 1 else g[(d, m, N, P)]
                    mm_, l_, h_, _ = stats(rs)
                    med.append(mm_); lo.append(l_); hi.append(h_)
                    S.append(T1 / mm_); E.append(T1 / mm_ / P)
                    share = statistics.median(r["t_comm"] / r["t_total"] for r in rs) if P > 1 else 0.0
                    line = f"| {d if P > 1 else 'serial'} | {m if P > 1 else '-'} | {P} | {fmt_t(mm_)} ({fmt_t(l_)}–{fmt_t(h_)}) | {T1 / mm_:.2f} | " \
                           f"{T1 / mm_ / P:.2f} | {share:.0%} |"
                    if model:
                        used = set(rs[0]["hosts_used"].split())
                        c, k = model.step(d if P > 1 else "serial", m, N, P, used)
                        c1, _ = model.step("serial", "-", N, 1, set(t1runs[0]["hosts_used"].split()))
                        line += f" {fmt_t(c + k)} | {c1 / (c + k):.2f} |"
                    if P > 1 or (d, m) == series[0]:
                        md.append(line)
                c = colors.get((d, m), None)
                err = [np.array(med) - np.array(lo), np.array(hi) - np.array(med)]
                axs[0].errorbar(Ps, med, yerr=err, fmt="o-", color=c, capsize=3, label=f"{d}, map-by {m}")
                axs[1].plot(Ps, S, "o-", color=c, label=f"{d}, map-by {m}")
                axs[2].plot(Ps, E, "o-", color=c, label=f"{d}, map-by {m}")
                if model:
                    c1, _ = model.step("serial", "-", N, 1, set(t1runs[0]["hosts_used"].split()))
                    used_by_P = {P: set(g[(d, m, N, P)][0]["hosts_used"].split()) for P in Ps if P > 1}
                    mt = [c1] + [sum(model.step(d, m, N, P, used_by_P[P])) for P in Ps[1:]]
                    axs[0].plot(Ps, mt, "--", color=c, alpha=0.8)
                    axs[1].plot(Ps, [c1 / t for t in mt], "--", color=c, alpha=0.8)
                    axs[2].plot(Ps, [c1 / t / P for t, P in zip(mt, Ps)], "--", color=c, alpha=0.8)
            md.append("")
            Pall = sorted({r["P"] for r in strong if r["N"] == N} | {1})
            axs[0].set_yscale("log")
            axs[0].set_ylabel("time per step [s]")
            axs[1].plot(Pall, Pall, ":", color="gray", label="ideal")
            axs[1].set_ylabel("speed-up S = T1 / TP")
            axs[2].axhline(1, ls=":", color="gray")
            axs[2].set_ylabel("efficiency E = S / P")
            for ax in axs:
                ax.set_xlabel("P")
                ax.set_xticks(Pall)
                ax.grid(True, alpha=0.3)
            axs[1].set_yscale("log", base=2)
            axs[0].legend(fontsize=8)
            fig.suptitle(f"Strong scaling, N = {N}" + ("   (markers: measured median with min–max; "
                                                         "dashed: model)" if model else ""))
            fig.tight_layout()
            fig.savefig(os.path.join(out, f"strong_N{N}.png"), dpi=150)
            plt.close(fig)

        # communication share plot
        fig, ax = plt.subplots(figsize=(7, 4.3))
        for (d, m, N), rs in sorted(group([r for r in strong if r["P"] > 1], ["decomp", "map", "N"]).items()):
            Ps = sorted({r["P"] for r in rs})
            sh = [statistics.median(r["t_comm"] / r["t_total"] for r in rs if r["P"] == P) for P in Ps]
            ax.plot(Ps, sh, "o-", label=f"{d}, {m}, N={N}")
        ax.set_xlabel("P")
        ax.set_xticks(sorted({r["P"] for r in strong if r["P"] > 1}))
        ax.set_ylabel("t_comm / t_total")
        ax.set_ylim(0, 1)
        ax.set_title("Share of time in the halo exchange (strong scaling)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "comm_share.png"), dpi=150)
        plt.close(fig)

        # crossover summary
        md += ["### Does spreading across VMs pay off?", "",
               "Best measured speed-up per N (any decomposition and placement, P > 1):", "",
               "| N | best S | at | model best S |", "|---|---|---|---|"]
        for N in sorted({r["N"] for r in strong}):
            t1 = g.get(("serial", "-", N, 1))
            if not t1:
                continue
            T1 = stats(t1)[0]
            best = max(((T1 / stats(rs)[0], k) for k, rs in g.items() if k[2] == N and k[3] > 1),
                       default=(float("nan"), None))
            mb = ""
            if model:
                c1, _ = model.step("serial", "-", N, 1, set(t1[0]["hosts_used"].split()))
                mb = max((c1 / sum(model.step(k[0], k[1], N, k[3], set(rs[0]["hosts_used"].split()))), k)
                         for k, rs in g.items() if k[2] == N and k[3] > 1)
                mb = f"{mb[0]:.2f} ({mb[1][0]}, {mb[1][1]}, P={mb[1][3]})"
            if best[1]:
                md.append(f"| {N} | {best[0]:.2f} | {best[1][0]}, map-by {best[1][1]}, P={best[1][3]} | {mb} |")
        md.append("")

    # ---- weak scaling
    weak = [r for r in rows if r["exp"] == "weak"]
    if weak:
        g = group(weak, ["decomp", "P"])
        md += ["## 3. Weak scaling", "",
               "Work per process fixed; ideally the time per step stays constant.", "",
               "| decomp | P | N | time/step (spread) | comm share |" + (" model |" if model else ""),
               "|---|---|---|---|---|" + ("---|" if model else "")]
        fig, ax = plt.subplots(figsize=(7, 4.3))
        serial = g.get(("serial", 1))
        allP = set()
        for d, col in (("strips", "C0"), ("blocks", "C2")):
            Ps = sorted(P for (dd, P) in g if dd == d)
            if not Ps:
                continue
            pts = ([(1, serial)] if serial else []) + [(P, g[(d, P)]) for P in Ps]
            xs, med, lo, hi, mt = [], [], [], [], []
            for P, rs in pts:
                m_, l_, h_, _ = stats(rs)
                xs.append(P); med.append(m_); lo.append(l_); hi.append(h_)
                share = statistics.median(r["t_comm"] / r["t_total"] for r in rs) if P > 1 else 0
                line = f"| {d if P > 1 else 'serial'} | {P} | {rs[0]['N']} | {fmt_t(m_)} ({fmt_t(l_)}–{fmt_t(h_)}) | {share:.0%} |"
                if model:
                    t = sum(model.step(d if P > 1 else "serial", "node", rs[0]["N"], P,
                                       set(rs[0]["hosts_used"].split())))
                    mt.append(t)
                    line += f" {fmt_t(t)} |"
                if P > 1 or d == "strips":
                    md.append(line)
            allP |= set(xs)
            ax.errorbar(xs, med, yerr=[np.array(med) - np.array(lo), np.array(hi) - np.array(med)],
                        fmt="o-", color=col, capsize=3, label=f"{d} (measured)")
            if model:
                ax.plot(xs, mt, "--", color=col, label=f"{d} (model)")
        md.append("")
        ax.set_xlabel("P")
        ax.set_ylabel("time per step [s]")
        ax.set_yscale("log")
        ax.set_xticks(sorted(allP))
        ax.set_title("Weak scaling (fixed work per process, map-by node)")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out, "weak.png"), dpi=150)
        plt.close(fig)

    # ---- real cluster prediction
    if model and strong:
        real = (a.real_ts, 1.0 / a.real_bw)
        Ns = sorted({r["N"] for r in strong})
        Pmax = max(r["P"] for r in strong)
        Ps = sorted({1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64} | set(range(1, Pmax + 1)))
        fig, ax = plt.subplots(figsize=(7, 4.3))
        md += ["## Prediction for a real cluster", "",
               f"Same compute rate, t_s = {fmt_t(real[0])}, bandwidth {a.real_bw / 1e9:.0f} GB/s "
               "(model only).", "", "| N | decomp | S at P=6 | S at P=64 |", "|---|---|---|---|"]
        for k, N in enumerate(Ns):
            c1 = model.step("serial", "-", N, 1)[0]
            for d, ls in (("strips", "o-"), ("blocks", "s--")):
                S = [c1 / sum(model.step(d, "node", N, P, None, real)) if P > 1 else 1.0 for P in Ps]
                ax.plot(Ps, S, ls, color=f"C{k}", ms=3, label=f"N={N}, {d}")
                md.append(f"| {N} | {d} | {S[Ps.index(6)]:.2f} | {S[Ps.index(64)]:.1f} |")
        md.append("")
        ax.plot(Ps, Ps, ":", color="gray", label="ideal")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.set_xlabel("P")
        ax.set_ylabel("model speed-up")
        ax.set_title(f"Model: real cluster (t_s = {fmt_t(real[0])}, {a.real_bw / 1e9:.0f} GB/s)")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "real_cluster.png"), dpi=150)
        plt.close(fig)

    failed = 0
    with open(os.path.join(a.folder, "runs.csv")) as f:
        failed = sum(1 for r in csv.DictReader(f) if r["ok"] != "1")
    head = [f"# Performance results: `{os.path.basename(os.path.normpath(a.folder))}`", "",
            f"{len(rows)} successful runs" + (f", {failed} failed (see runs.csv)" if failed else "") + ".", ""]
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(head + md) + "\n")
    print("\n".join(head + md))
    print(f"\nwrote {out}/")


if __name__ == "__main__":
    main()
