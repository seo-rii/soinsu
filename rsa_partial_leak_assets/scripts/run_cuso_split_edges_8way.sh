#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/ct07_cuso_$(date +%Y%m%d_%H%M%S)/split_edges}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_GRAPH="${CUSO_GRAPH:-on}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-0}"
CUSO_UPPER_LOW_VARIABLES="${CUSO_UPPER_LOW_VARIABLES:-0}"
CUSO_UPPER_ALL_VARIABLES="${CUSO_UPPER_ALL_VARIABLES:-0}"
CUSO_UPPER_VARIABLES="${CUSO_UPPER_VARIABLES:-}"
mkdir -p "$LOG_DIR"

origin_args=()
if [[ "$CUSO_UPPER_LOW_VARIABLES" == "1" || "$CUSO_UPPER_LOW_VARIABLES" == "true" ]]; then
  origin_args+=(--cuso-upper-low-variables)
fi
if [[ "$CUSO_UPPER_ALL_VARIABLES" == "1" || "$CUSO_UPPER_ALL_VARIABLES" == "true" ]]; then
  origin_args+=(--cuso-upper-all-variables)
fi
for spec in $CUSO_UPPER_VARIABLES; do
  origin_args+=(--cuso-upper-variable "$spec")
done

pids=()
for i in $(seq 0 7); do
  A=$((i*32))
  B=$(((i+1)*32))
  (
    log="$LOG_DIR/split_${A}_${B}.log"
    echo "edge_range=$A:$B graph=$CUSO_GRAPH upper_low=$CUSO_UPPER_LOW_VARIABLES upper_all=$CUSO_UPPER_ALL_VARIABLES upper_variables=$CUSO_UPPER_VARIABLES" | tee "$log"
    set +e
    cmd=(sage -python src/solve7_main.py \
      --mode cuso-split \
      --cuso-split-brute-small-edges \
      --a "$A" --b "$B" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph "$CUSO_GRAPH" \
      "${origin_args[@]}")
    if [[ "$CUSO_TIMEOUT" == "0" || -z "$CUSO_TIMEOUT" ]]; then
      "${cmd[@]}" 2>&1 | tee -a "$log"
    else
      timeout "$CUSO_TIMEOUT" "${cmd[@]}" 2>&1 | tee -a "$log"
    fi
    status=${PIPESTATUS[0]}
    echo "exit_status=$status" | tee -a "$log"
  ) &
  pids+=("$!")
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

grep -RniE "FOUND|p =|q =|m.hex|split roots|Traceback|Exception|cuso unavailable|cuso fail|not found|exit_status" "$LOG_DIR" || true
echo "logs written to $LOG_DIR"
