#!/usr/bin/env python3
"""convergence.py -- convergence study for Task 7.2 of the background notes.

Runs the verification case V to T = 0.05 for N = 16, 32, 64, 128, 256 with r
fixed, collects ||e||_inf and ||e||_2h from the solver's summary line, computes
the observed order p_obs = log2(E(h) / E(h/2)) and writes:

    <outdir>/convergence.md    the table for the report
    <outdir>/convergence.png   log-log plot of the errors against h
    <outdir>/raw.txt           the solver's raw output lines

Usage (from the repository root):
    python3 scripts/convergence.py [--r 0.2] [--T 0.05] [--outdir results/convergence]

r = 0.2 is the default because T/dt = 0.05 N^2 / r = N^2 / 4 is then an
integer for every even N, so r is exactly the same on every grid.
"""
import argparse
import math
import os
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEAT = os.path.join(ROOT, "build", "heat")


def run(N, r, T):
    cmd = [HEAT, "--case", "V", "--N", str(N), "--r", str(r), "--T", str(T), "--decomp", "serial"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
    kv = dict(item.split("=", 1) for item in out.split())
    return out, kv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", type=float, default=0.2)
    ap.add_argument("--T", type=float, default=0.05)
    ap.add_argument("--Ns", type=int, nargs="+", default=[16, 32, 64, 128, 256])
    ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "convergence"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows, raw = [], []
    for N in a.Ns:
        line, kv = run(N, a.r, a.T)
        raw.append(line)
        print(line)
        rows.append(dict(N=N, h=1.0 / N, steps=int(kv["steps"]), r=float(kv["r"]),
                         einf=float(kv["err_inf"]), e2=float(kv["err_2h"])))

    lines = ["| N | h | steps | r used | ‖e‖∞ | ‖e‖₂,h | p_obs (∞) | p_obs (2,h) |",
             "|---|---|---|---|---|---|---|---|"]
    for k, row in enumerate(rows):
        if k == 0:
            p_inf = p_2 = "—"
        else:
            p_inf = f"{math.log2(rows[k-1]['einf'] / row['einf']):.3f}"
            p_2 = f"{math.log2(rows[k-1]['e2'] / row['e2']):.3f}"
        lines.append(f"| {row['N']} | {row['h']:.6g} | {row['steps']} | {row['r']:.6g} | "
                     f"{row['einf']:.4e} | {row['e2']:.4e} | {p_inf} | {p_2} |")
    table = "\n".join(lines)
    print()
    print(table)

    with open(os.path.join(a.outdir, "convergence.md"), "w") as f:
        f.write(f"Verification case V, T = {a.T}, r = {a.r}\n\n{table}\n")
    with open(os.path.join(a.outdir, "raw.txt"), "w") as f:
        f.write("\n".join(raw) + "\n")

    h = [row["h"] for row in rows]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(h, [row["e2"] for row in rows], "o-", label=r"$\|e\|_{2,h}$")
    ax.loglog(h, [row["einf"] for row in rows], "s-", label=r"$\|e\|_\infty$")
    ref = [rows[0]["e2"] * (x / h[0]) ** 2 for x in h]
    ax.loglog(h, ref, "k--", lw=1, label=r"slope 2 reference")
    ax.set_xlabel("h")
    ax.set_ylabel("error at T")
    ax.set_title(f"Convergence study, case V, T={a.T}, r={a.r}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "convergence.png"), dpi=150)
    print(f"\nwrote {a.outdir}/convergence.md, convergence.png, raw.txt")


if __name__ == "__main__":
    main()
