#!/usr/bin/env bash
# run_pingpong.sh -- run the ping-pong benchmark for every VM on its own
# (two processes on the same VM) and for every pair of VMs, and record
# whether each pair talks directly or through a Tailscale DERP relay.
#
# Usage (from the repository root, on any one VM):
#   scripts/run_pingpong.sh heat-lynn heat-ghadi heat-rim
#
# Extra arguments for the pingpong program can be passed via PP_ARGS, e.g.
#   PP_ARGS="--max 4194304" scripts/run_pingpong.sh heat-lynn heat-ghadi
#
# Requirements: build/pingpong exists at the SAME path on every VM
# (clone the repo to the same place and run `make` on each VM), and
# passwordless SSH works between all VMs.
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "usage: $0 host1 [host2 ...]" >&2
    exit 1
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/build/pingpong"
OUT="$ROOT/results/pingpong/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT"

NET=(--mca oob_tcp_if_include tailscale0 --mca btl_tcp_if_include tailscale0)
# Same-VM runs use shared memory (vader) by default, which is how two
# processes on one machine normally talk. Set SAME_VM_BTL=tcp,self to force
# TCP instead and compare.
SAME_VM_BTL="${SAME_VM_BTL:-self,vader}"
PP_ARGS="${PP_ARGS:-}"

hosts=("$@")

{
    echo "date: $(date -Is)"
    echo "launched from: $(hostname)"
    echo "hosts: ${hosts[*]}"
    echo "same-VM btl: $SAME_VM_BTL"
    echo "pingpong args: $PP_ARGS"
    mpirun --version | head -1
} > "$OUT/info.txt"

echo "Results go to $OUT"

for h in "${hosts[@]}"; do
    echo "== same VM: $h"
    # shellcheck disable=SC2086
    mpirun -np 2 --host "$h:2" "${NET[@]}" --mca btl "$SAME_VM_BTL" \
        "$BIN" $PP_ARGS | tee "$OUT/${h}__${h}.csv"
done

for ((i = 0; i < ${#hosts[@]}; i++)); do
    for ((j = i + 1; j < ${#hosts[@]}; j++)); do
        a="${hosts[i]}"; b="${hosts[j]}"
        echo "== path check: $a -> $b"
        ssh "$a" tailscale ping -c 5 "$b" 2>&1 | tee "$OUT/${a}__${b}.path.txt" || true
        echo "== cross VM: $a <-> $b"
        # shellcheck disable=SC2086
        mpirun -np 2 --host "$a,$b" "${NET[@]}" --mca btl tcp,self \
            "$BIN" $PP_ARGS | tee "$OUT/${a}__${b}.csv"
    done
done

echo
echo "Done. Fit the model with:"
echo "  python3 scripts/fit_pingpong.py $OUT"
