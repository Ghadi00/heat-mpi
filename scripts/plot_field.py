#!/usr/bin/env python3
"""plot_field.py -- draw a heat map of a field file written by `heat --out`.

Usage:
    python3 scripts/plot_field.py u.bin [--title "Case 2, N=256"] [--out u.png]

The file holds (N+1)^2 little-endian doubles, row-major, row j = line y = y_j;
N is worked out from the file size.
"""
import argparse
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--symmetric", action="store_true",
                    help="colour scale symmetric around 0 (useful for unstable runs)")
    a = ap.parse_args()

    u = np.fromfile(a.field, dtype="<f8")
    n1 = math.isqrt(u.size)
    if n1 * n1 != u.size:
        raise SystemExit(f"{a.field}: {u.size} values is not a square grid")
    u = u.reshape(n1, n1)          # u[j, i]: row j is y = y_j
    N = n1 - 1

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    kw = {}
    if a.symmetric:
        m = np.nanmax(np.abs(u))
        kw = dict(cmap="RdBu_r", vmin=-m, vmax=m)
    else:
        kw = dict(cmap="inferno")
    im = ax.imshow(u, origin="lower", extent=[0, 1, 0, 1], **kw)
    fig.colorbar(im, ax=ax, label="u")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(a.title or f"{os.path.basename(a.field)}  (N={N})")
    fig.tight_layout()
    out = a.out or os.path.splitext(a.field)[0] + ".png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}  (min={u.min():.4g}, max={u.max():.4g})")


if __name__ == "__main__":
    main()
