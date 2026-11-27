# Parallel 2D heat diffusion with MPI over Tailscale

## Build

On every VM, clone this repository to the **same path** (e.g. `~/heat-mpi`)
and run:

    make

This builds `build/heat` (the solver), `build/hello` and `build/pingpong`.
Flags: `-O3 -std=c++17 -ffp-contract=off` (no fused multiply-adds, so results
are bit-for-bit reproducible).

## The solver: `build/heat`

    build/heat --case 2 --N 1024 --r 0.2 --steps 500 --out u.bin

| Option | Meaning | Example |
|---|---|---|
| `--case` | `1` (four hot edges), `2` (hot top and bottom), `V` (verification) | `--case 2` |
| `--N` | number of intervals per direction (grid has (N+1)² nodes) | `--N 1024` |
| `--r` | α∆t/h². Must be ≤ 1/4; larger values are refused | `--r 0.2` |
| `--T` | final time; ∆t is chosen so T/∆t is an integer and r is not exceeded | `--T 0.05` |
| `--steps` | fixed number of steps instead of `--T` (∆t = r h²) | `--steps 500` |
| `--out` | write the final field: (N+1)² raw little-endian doubles, row-major, row j = line y = y_j, no header | `--out u.bin` |
| `--decomp` | `serial` (1 process, reference), `strips` (Part 2, 1D row strips) or `blocks` (Part 3, Px × Py grid chosen by `MPI_Dims_create`; printed as `Px=`/`Py=` in the summary) | `--decomp blocks` |
| `--halo` | strips only (blocks always use the shift order): `shift` (default; every rank sends up, then every rank sends down) or `pair` (exchange with the upper neighbour, then the lower; the order written in the brief, slower for large P) | `--halo pair` |
| `--print` | debug: print every rank's local array including ghost rows, in rank order | |
| `--allow-unstable` | run even if r > 1/4 (only for Task 6.7) | |

Give exactly one of `--T` and `--steps`. Output is one line of `key=value`
pairs, e.g.

    decomp=serial P=1 case=V N=64 steps=1024 r=0.2 dt=4.882813e-05 T=0.05 t_total=0.002 t_comp=0.002 t_comm=0 err_inf=1.03e-04 err_2h=5.17e-05

(`err_inf`, `err_2h` only for case V). Times cover the time-stepping loop only;
in parallel each is the maximum over ranks. With `--decomp strips` or `blocks`, the
field is gathered on rank 0 only when `--out` is given.

Running the strip version across the cluster:

    mpirun --hostfile hosts -np 6 --map-by node \
        --mca oob_tcp_if_include tailscale0 --mca btl_tcp_if_include tailscale0 \
        --mca btl tcp,self \
        build/heat --case 2 --N 512 --r 0.24 --steps 100 --decomp strips

Helpers:

    python3 scripts/convergence.py            # Task 7.2 table + plot -> results/convergence/
    scripts/check_correctness.sh              # serial vs strips and blocks, bit-for-bit (cmp), on this VM
    CLUSTER=1 scripts/check_correctness.sh    # the same across the VMs in ./hosts
    python3 scripts/plot_field.py u.bin       # heat map -> u.png

## The cluster

`hosts.example` lists the VMs and their slots (= vCPU count). Copy it to
`hosts` and adjust. Every mpirun across VMs uses the Tailscale interface:

    mpirun --hostfile hosts -np 6 \
        --mca oob_tcp_if_include tailscale0 \
        --mca btl_tcp_if_include tailscale0 \
        --mca btl tcp,self \
        build/hello

## Ping-pong (Part 0, Step 6)

    scripts/run_pingpong.sh heat-lynn heat-ghadi heat-rim
    python3 scripts/fit_pingpong.py results/pingpong/<run-folder>

`build/pingpong` options:

| Option | Meaning | Default |
|---|---|---|
| `--max BYTES` | largest message size | 16777216 (16 MiB) |
| `--budget SEC` | measuring time per message size | 2.0 |
| `--min-reps N` | minimum repetitions per size | 5 |
| `--max-reps N` | maximum repetitions per size | 1000 |

## LLM use

See `CLAUDE.md` (also used as `AGENTS.md`) and the `llm_log/` folder.

## Performance study (Part 4)

Run on the VM listed first in `hosts`, after `git pull && make` on every VM:

    python3 scripts/perf_run.py --dry-run --pingpong results/pingpong/<folder>   # plan + time estimate
    python3 scripts/perf_run.py --tag <network-name>                             # run everything
    python3 scripts/perf_analyze.py results/perf/<folder> --pingpong results/pingpong/<folder>

`perf_run.py` runs the three experiments of Section 8.2 of the brief (baseline,
strong scaling with `--map-by node` and `--map-by slot`, weak scaling), each
configuration `--reps` times (default 3), and appends every result immediately to
`results/perf/<folder>/runs.csv` together with the exact command. A stopped session
continues with `--resume <folder>`; a hung run is killed after `--timeout` seconds;
after `--max-fail` failures in a row (default 5) it stops. `network.log` records
`tailscale status` and the path of every pair at the start of each experiment.
Write the network setup by hand into `<folder>/NETWORK.md`.

| Option | Default | Meaning |
|---|---|---|
| `--exp` | `baseline strong weak` | which experiments |
| `--steps` | 100 | steps for strong and weak scaling |
| `--baseline-steps` | 500 | steps for the serial baseline |
| `--reps` | 3 | repetitions of every configuration |
| `--Ns` | `512 2048 8192` | strong-scaling grid sizes |
| `--Ps` | 1 .. all slots | process counts |
| `--weak-per` | 1024 | interior nodes per process = this² |
| `--timeout` | 1800 | seconds before a run is killed |

`perf_analyze.py` takes the median and spread of the repetitions and writes
`analysis/summary.md` (all tables) and the plots: compute rate per VM, strong
scaling (time per step, speed-up, efficiency, measured and model), weak scaling,
communication share, and the model's prediction for a real cluster. The model is
described at the top of the script.
