#!/usr/bin/env bash
# check_correctness.sh -- run the serial solver and the parallel
# decompositions for several N and P and check with `cmp` that every parallel
# output file is bit-for-bit identical to the serial one (deliverable 3).
#
# Usage (from the repository root):
#   scripts/check_correctness.sh                       # on this VM only
#   CLUSTER=1 scripts/check_correctness.sh             # across the VMs in ./hosts
#
# Settings (environment variables, defaults in brackets):
#   CASES    test cases                    ["1 2 V"]
#   NS       grid sizes N                  ["8 17 64 101"]
#            (N-1 = 7, 16, 63, 100: divisible and not divisible by P = 2, 3, 4)
#   PS       process counts                ["1 2 3 4"]
#   DECOMPS  parallel decompositions       ["strips"]   (add "blocks" in Part 3)
#   HALOS    halo variants for strips      ["shift pair"]
#   STEPS    time steps per run            [50]
#   R        alpha dt / h^2                [0.24]
#   CLUSTER  1 = use ./hosts and Tailscale, spreading ranks over the VMs
#            (--map-by node); default: all ranks on this VM (--oversubscribe)
#   EXTRA    extra mpirun arguments        [""]
#
# Output: results/correctness/<date-time>/log.txt and summary.txt.
# Exit status 0 only if every comparison passed.
set -uo pipefail
trap 'echo; echo "interrupted"; pkill -P $$ mpirun 2>/dev/null; exit 130' INT TERM

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HEAT="$ROOT/build/heat"
CASES="${CASES:-1 2 V}"
NS="${NS:-8 17 64 101}"
PS="${PS:-1 2 3 4}"
DECOMPS="${DECOMPS:-strips}"
HALOS="${HALOS:-shift pair}"
STEPS="${STEPS:-50}"
R="${R:-0.24}"
EXTRA="${EXTRA:-}"

if [ "${CLUSTER:-0}" = "1" ]; then
    MPIRUN=(mpirun --hostfile "$ROOT/hosts" --map-by node
            --mca oob_tcp_if_include tailscale0 --mca btl_tcp_if_include tailscale0
            --mca btl tcp,self)
    WHERE="cluster ($(grep -v '^#' "$ROOT/hosts" | awk '{print $1}' | xargs))"
else
    MPIRUN=(mpirun --oversubscribe)
    WHERE="this VM only ($(hostname))"
fi

OUT="$ROOT/results/correctness/$(date +%Y%m%d-%H%M%S)"
TMP="$OUT/fields"
mkdir -p "$TMP"
LOG="$OUT/log.txt"
SUM="$OUT/summary.txt"

{
    echo "date: $(date -Is)"
    echo "where: $WHERE"
    echo "commit: $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    echo "cases=[$CASES] N=[$NS] P=[$PS] decomps=[$DECOMPS] halos=[$HALOS] steps=$STEPS r=$R"
    echo
} | tee "$SUM" > "$LOG"

pass=0; fail=0
printf "%-5s %-5s %-8s %-6s %-3s %s\n" case N decomp halo P result | tee -a "$SUM"

for c in $CASES; do
    for N in $NS; do
        ref="$TMP/serial_c${c}_N${N}.bin"
        args=(--case "$c" --N "$N" --r "$R" --steps "$STEPS")
        echo "\$ heat ${args[*]} --decomp serial" >> "$LOG"
        if ! "$HEAT" "${args[@]}" --decomp serial --out "$ref" >> "$LOG" 2>&1; then
            echo "serial run failed for case $c N=$N" | tee -a "$SUM"; fail=$((fail + 1)); continue
        fi
        for d in $DECOMPS; do
            halos="-"; [ "$d" = "strips" ] && halos="$HALOS"
            for halo in $halos; do
                for P in $PS; do
                    if [ "$d" = "strips" ] && [ "$P" -gt $((N - 1)) ]; then continue; fi
                    f="$TMP/${d}_${halo}_c${c}_N${N}_P${P}.bin"
                    hargs=(); [ "$halo" != "-" ] && hargs=(--halo "$halo")
                    echo "\$ mpirun -np $P heat ${args[*]} --decomp $d ${hargs[*]}" >> "$LOG"
                    # shellcheck disable=SC2086
                    if "${MPIRUN[@]}" $EXTRA -np "$P" "$HEAT" "${args[@]}" --decomp "$d" \
                           "${hargs[@]}" --out "$f" >> "$LOG" 2>&1 && cmp -s "$ref" "$f"; then
                        res="IDENTICAL"; pass=$((pass + 1))
                    else
                        res="DIFFERENT"; fail=$((fail + 1))
                        cmp "$ref" "$f" >> "$LOG" 2>&1
                    fi
                    printf "%-5s %-5s %-8s %-6s %-3s %s\n" "$c" "$N" "$d" "$halo" "$P" "$res" | tee -a "$SUM"
                done
            done
        done
    done
done

echo | tee -a "$SUM"
echo "passed: $pass   failed: $fail" | tee -a "$SUM"
rm -rf "$TMP"          # the .bin files are large and reproducible; keep the logs
echo "logs in $OUT"
[ "$fail" -eq 0 ]
