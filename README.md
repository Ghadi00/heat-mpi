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
| `--decomp` | `serial` (strips and blocks: Parts 2 and 3) | `--decomp serial` |
| `--allow-unstable` | run even if r > 1/4 (only for Task 6.7) | |

Give exactly one of `--T` and `--steps`. Output is one line of `key=value`
pairs, e.g.

    decomp=serial P=1 case=V N=64 steps=1024 r=0.2 dt=4.882813e-05 T=0.05 t_total=0.002 t_comp=0.002 t_comm=0 err_inf=1.03e-04 err_2h=5.17e-05

(`err_inf`, `err_2h` only for case V). Times cover the time-stepping loop only.

Helpers:

    python3 scripts/convergence.py            # Task 7.2 table + plot -> results/convergence/
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
