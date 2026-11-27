#!/usr/bin/env python3
"""perf_run.py -- run the Part 4 performance experiments and keep every raw result.

Experiments (project brief, Section 8.2):
  baseline  serial run on EACH VM for N = 512 1024 2048 4096 8192
            -> compute rate (time per node per step) of every machine
  strong    N = 512 2048 8192, strips and blocks, P = 1 .. all slots,
            placements --map-by node (round-robin) and --map-by slot (fill a VM first)
            P = 1 is the serial version on the first host (the reference T1)
  weak      fixed work per process (default 1024 x 1024 interior nodes), N grows with P,
            strips and blocks, --map-by node
Communication share is not a separate experiment: every run records t_comm / t_total.

Usage (from the repository root, on the VM listed first in ./hosts):
    python3 scripts/perf_run.py --dry-run                 # list what would run, nothing runs
    python3 scripts/perf_run.py                           # everything, default settings
    python3 scripts/perf_run.py --exp baseline            # only one experiment
    python3 scripts/perf_run.py --exp strong --Ns 512 2048 --steps 200
    python3 scripts/perf_run.py --resume results/perf/<folder>   # carry on after a stop

Every finished run is appended at once to <folder>/runs.csv, so nothing is lost if the
script is stopped; --resume skips runs already in runs.csv (failed ones are retried).
A run that takes longer than --timeout seconds is killed and recorded as failed.

Rules from the brief that this script enforces:
  * fixed number of steps per experiment, r close to the stability limit (0.24)
  * every configuration repeated --reps times (default 3); the median is taken later
  * never more processes on a VM than its slots (= vCPUs) in ./hosts (no --oversubscribe)
  * raw output kept (runs.csv has every summary line and the exact command)

Analyse with:  python3 scripts/perf_analyze.py results/perf/<folder>
"""
import argparse
import csv
import datetime as dt
import math
import os
import shlex
import subprocess
import sys
import time
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEAT = os.path.join(ROOT, "build", "heat")
FIELDS = ["exp", "decomp", "map", "N", "P", "rep", "steps", "ok",
          "t_total", "t_comp", "t_comm", "Px", "Py", "host", "hosts_used",
          "started", "wall_s", "cmd", "output"]


def read_hosts(path):
    hosts = []
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            slots = 1
            for p in parts[1:]:
                if p.startswith("slots="):
                    slots = int(p.split("=")[1])
            hosts.append((parts[0], slots))
    if not hosts:
        sys.exit(f"no hosts in {path}")
    return hosts


def hosts_used(hosts, P, mapping):
    """Which VMs Open MPI will use for P ranks with --map-by node / slot."""
    used = []
    if mapping == "slot":
        left = P
        for h, s in hosts:
            if left <= 0:
                break
            used.append(h)
            left -= s
    else:
        used = [h for h, _ in hosts[:P]]
    return used


def plan(a, hosts):
    runs = []   # dicts with exp, decomp, map, N, P, steps, host (for baseline)
    total_slots = sum(s for _, s in hosts)
    Pmax = a.Pmax or total_slots
    Ps = [p for p in (a.Ps or range(1, Pmax + 1)) if p <= total_slots]
    if "baseline" in a.exp:
        for N in a.baseline_Ns:
            for h, _ in hosts:
                runs.append(dict(exp="baseline", decomp="serial", map="-", N=N, P=1,
                                 steps=a.baseline_steps, host=h))
    if "strong" in a.exp:
        for N in a.Ns:
            runs.append(dict(exp="strong", decomp="serial", map="-", N=N, P=1,
                             steps=a.steps, host=hosts[0][0]))
            for d in a.decomps:
                for m in a.maps:
                    for P in Ps:
                        if P == 1:
                            continue
                        runs.append(dict(exp="strong", decomp=d, map=m, N=N, P=P,
                                         steps=a.steps, host=""))
    if "weak" in a.exp:
        for P in Ps:
            # (N-1)^2 interior nodes = P * per^2  ->  N = per*sqrt(P) + 1, rounded
            N = int(round(a.weak_per * math.sqrt(P))) + 1
            if P == 1:
                runs.append(dict(exp="weak", decomp="serial", map="-", N=N, P=1,
                                 steps=a.weak_steps or a.steps, host=hosts[0][0]))
                continue
            for d in a.decomps:
                runs.append(dict(exp="weak", decomp=d, map="node", N=N, P=P,
                                 steps=a.weak_steps or a.steps, host=""))
    return runs


def command(run, a, hostfile):
    NET = [] if a.iface == "none" else ["--mca", "oob_tcp_if_include", a.iface,
                                        "--mca", "btl_tcp_if_include", a.iface]
    heat = [HEAT, "--case", "2", "--N", str(run["N"]), "--r", str(a.r),
            "--steps", str(run["steps"]), "--decomp", run["decomp"]]
    if run["decomp"] == "serial":
        return ["mpirun", "-np", "1", "--host", run["host"]] + NET + heat
    return (["mpirun", "--hostfile", hostfile, "-np", str(run["P"]), "--map-by", run["map"]]
            + NET + ["--mca", "btl", "tcp,self,vader"] + heat)


def key(r):
    return (r["exp"], r["decomp"], r["map"], str(r["N"]), str(r["P"]), str(r["steps"]),
            r.get("host", ""), str(r["rep"]))


def parse_summary(text):
    for line in text.splitlines():
        if line.startswith("decomp="):
            return dict(kv.split("=", 1) for kv in line.split() if "=" in kv), line
    return None, ""


def snapshot(folder, label, hosts):
    """Record the network state (tailscale status and paths) into network.log."""
    with open(os.path.join(folder, "network.log"), "a") as f:
        f.write(f"\n===== {label}  {dt.datetime.now().isoformat(timespec='seconds')}\n")
        try:
            f.write(subprocess.run(["tailscale", "status"], capture_output=True, text=True,
                                   timeout=20).stdout)
        except Exception as e:     # tailscale missing (local test) or hanging
            f.write(f"(tailscale status failed: {e})\n")
        names = [h for h, _ in hosts]
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                try:
                    out = subprocess.run(["ssh", names[i], "tailscale", "ping", "-c", "2", names[j]],
                                         capture_output=True, text=True, timeout=30)
                    f.write(f"-- {names[i]} -> {names[j]}\n{out.stdout}{out.stderr}")
                except Exception as e:
                    f.write(f"-- {names[i]} -> {names[j]}: {e}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exp", nargs="+", default=["baseline", "strong", "weak"],
                    choices=["baseline", "strong", "weak"])
    ap.add_argument("--hosts", default=os.path.join(ROOT, "hosts"))
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--steps", type=int, default=100,
                    help="steps for strong and weak scaling (default 100: over a slow network one "
                         "N=8192 step can take ~1 s; time per step is what is compared)")
    ap.add_argument("--baseline-steps", type=int, default=500,
                    help="steps for the serial baseline (no network involved, so it can be longer)")
    ap.add_argument("--weak-steps", type=int, default=None)
    ap.add_argument("--r", type=float, default=0.24)
    ap.add_argument("--Ns", type=int, nargs="+", default=[512, 2048, 8192], help="strong scaling grid sizes")
    ap.add_argument("--baseline-Ns", type=int, nargs="+", default=[512, 1024, 2048, 4096, 8192])
    ap.add_argument("--weak-per", type=int, default=1024, help="interior nodes per process = this^2")
    ap.add_argument("--Ps", type=int, nargs="+", default=None, help="process counts (default 1..all slots)")
    ap.add_argument("--Pmax", type=int, default=None)
    ap.add_argument("--decomps", nargs="+", default=["strips", "blocks"])
    ap.add_argument("--maps", nargs="+", default=["node", "slot"])
    ap.add_argument("--timeout", type=float, default=1800, help="seconds before one run is killed")
    ap.add_argument("--resume", default=None, help="existing results folder to continue")
    ap.add_argument("--tag", default="", help="text added to the folder name, e.g. home-wifi")
    ap.add_argument("--iface", default="tailscale0", help="network interface for MPI ('none' = let MPI choose)")
    ap.add_argument("--max-fail", type=int, default=5,
                    help="stop after this many failed runs in a row (a dead link or a sleeping VM)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--pingpong", default=None,
                    help="with --dry-run: ping-pong folder used to ESTIMATE the total run time")
    ap.add_argument("--tnode", type=float, default=2e-9,
                    help="with --dry-run: guessed serial time per node per step [s] for the estimate")
    ap.add_argument("--no-network-log", action="store_true", help="skip tailscale snapshots (local tests)")
    a = ap.parse_args()

    hosts = read_hosts(a.hosts)
    runs = plan(a, hosts)
    todo = [dict(r, rep=k) for k in range(1, a.reps + 1) for r in runs]

    if a.dry_run:
        from collections import Counter
        c = Counter((r["exp"], r["decomp"], r["map"]) for r in todo)
        print(f"hosts: {hosts}")
        for (e, d, m), n in sorted(c.items()):
            print(f"  {e:8s} {d:7s} map={m:5s} {n:4d} runs")
        print(f"total {len(todo)} runs ({len(runs)} configurations x {a.reps} repetitions)")
        if a.pingpong:
            from perf_analyze import Model, pingpong_params
            net, loc, _, _ = pingpong_params(a.pingpong, "worst")
            mdl = Model({}, net, loc, hosts[0][1])
            mdl.t_node = lambda N, h: a.tnode
            per = defaultdict(float)
            for r in todo:
                used = [r["host"]] if r["decomp"] == "serial" else hosts_used(hosts, r["P"], r["map"])
                t = sum(mdl.step(r["decomp"], r["map"], r["N"], r["P"], set(used))) * r["steps"]
                per[r["exp"]] += t + 3.0           # + ~3 s to start mpirun
            print(f"\nESTIMATE (t_s={net[0]*1e3:.0f} ms, {1/net[1]/1e6:.2f} MB/s, "
                  f"t_node={a.tnode*1e9:.1f} ns -- rough!):")
            for e, t in per.items():
                print(f"  {e:8s} ~{t / 3600:.1f} h")
            print(f"  total    ~{sum(per.values()) / 3600:.1f} h")
        for r in todo[:: max(1, len(todo) // 8)]:
            print("  e.g.", " ".join(shlex.quote(x) for x in command(r, a, a.hosts)))
        return

    if not os.path.exists(HEAT):
        sys.exit(f"{HEAT} not found: run make first")

    if a.resume:
        folder = a.resume
    else:
        stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        folder = os.path.join(ROOT, "results", "perf", stamp + (f"-{a.tag}" if a.tag else ""))
    os.makedirs(folder, exist_ok=True)
    csv_path = os.path.join(folder, "runs.csv")

    done = set()
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                if row["ok"] == "1":
                    done.add(key(row))
    new_file = not os.path.exists(csv_path)

    with open(os.path.join(folder, "info.txt"), "a") as f:
        f.write(f"started {dt.datetime.now().isoformat(timespec='seconds')} on {os.uname().nodename}\n")
        f.write("args: " + " ".join(shlex.quote(x) for x in sys.argv[1:]) + "\n")
        try:
            commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
                                    capture_output=True, text=True).stdout.strip()
        except Exception:
            commit = "unknown"
        f.write(f"commit {commit}\nhosts {hosts}\n")
        f.write("Write the network setup (who was on which network, Wi-Fi or wired) in NETWORK.md.\n")

    if not a.no_network_log:
        snapshot(folder, "start", hosts)

    todo = [r for r in todo if key(dict(r, N=str(r["N"]), P=str(r["P"]), steps=str(r["steps"]),
                                        rep=str(r["rep"]))) not in done]
    print(f"results: {folder}")
    print(f"{len(todo)} runs to do ({len(done)} already done)")

    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()
        last_exp = None
        fails_in_row = 0
        for n, r in enumerate(todo, 1):
            if r["exp"] != last_exp and not a.no_network_log and last_exp is not None:
                snapshot(folder, f"before {r['exp']}", hosts)
            last_exp = r["exp"]
            cmd = command(r, a, a.hosts)
            started = dt.datetime.now().isoformat(timespec="seconds")
            t0 = time.time()
            try:
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout)
                text = p.stdout + p.stderr
                ok = p.returncode == 0
            except KeyboardInterrupt:
                subprocess.run(["pkill", "-u", os.environ.get("USER", ""), "-x", "mpirun"])
                print(f"\ninterrupted: the {n - 1} finished runs are saved. Continue later with\n"
                      f"  python3 scripts/perf_run.py --resume {folder} <same options as before>", flush=True)
                sys.exit(130)
            except subprocess.TimeoutExpired:
                text = f"TIMEOUT after {a.timeout:.0f} s"
                ok = False
                subprocess.run(["pkill", "-u", os.environ.get("USER", ""), "-x", "mpirun"])
            wall = time.time() - t0
            s, line = parse_summary(text)
            ok = ok and s is not None
            row = dict(exp=r["exp"], decomp=r["decomp"], map=r["map"], N=r["N"], P=r["P"],
                       rep=r["rep"], steps=r["steps"], ok=int(ok),
                       t_total=s.get("t_total", "") if s else "", t_comp=s.get("t_comp", "") if s else "",
                       t_comm=s.get("t_comm", "") if s else "", Px=s.get("Px", "") if s else "",
                       Py=s.get("Py", "") if s else "", host=r.get("host", ""),
                       hosts_used=" ".join([r["host"]] if r["decomp"] == "serial"
                                           else hosts_used(hosts, r["P"], r["map"])),
                       started=started, wall_s=f"{wall:.1f}",
                       cmd=" ".join(shlex.quote(x) for x in cmd),
                       output=line if ok else text.strip().replace("\n", " | ")[:2000])
            w.writerow(row)
            f.flush()
            status = (f"t_total={row['t_total']}s comm={float(row['t_comm']) / max(float(row['t_total']), 1e-12):.0%}"
                      if ok else "FAILED")
            where = r.get("host") or row["hosts_used"]
            fails_in_row = 0 if ok else fails_in_row + 1
            print(f"[{n}/{len(todo)}] {r['exp']:8s} {r['decomp']:7s} map={r['map']:5s} "
                  f"N={r['N']:5d} P={r['P']} rep={r['rep']} on {where}: {status}  ({wall:.0f}s wall)",
                  flush=True)
            if not ok:
                print("    " + row["output"][:300], flush=True)
            if fails_in_row >= a.max_fail:
                print(f"\n{fails_in_row} runs failed in a row: stopping. Check that every VM is awake "
                      f"(tailscale status, ssh <host> hostname), then continue with\n"
                      f"  python3 scripts/perf_run.py --resume {folder} <same options as before>")
                sys.exit(1)

    if not a.no_network_log:
        snapshot(folder, "end", hosts)
    print(f"\ndone. results in {folder}\nanalyse: python3 scripts/perf_analyze.py {folder}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\ninterrupted -- finished runs are saved; continue with --resume <folder>")
        sys.exit(130)
