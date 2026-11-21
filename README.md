# Parallel 2D heat diffusion with MPI over Tailscale

## Build

On every VM, clone this repository to the **same path** (e.g. `~/heat-mpi`)
and run:

    make

This builds `build/hello` and `build/pingpong` (the solver `build/heat`
will be added in Part 1).

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
