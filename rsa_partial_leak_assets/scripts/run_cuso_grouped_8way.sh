#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

LOG_DIR="${1:-logs/ct07_cuso_$(date +%Y%m%d_%H%M%S)/grouped}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_GRAPH="${CUSO_GRAPH:-on}"
mkdir -p "$LOG_DIR"

pids=()
for i in $(seq 0 7); do
  A=$((i*32))
  B=$(((i+1)*32))
  (
    log="$LOG_DIR/grouped_${A}_${B}.log"
    echo "range=$A:$B graph=$CUSO_GRAPH" | tee "$log"
    set +e
    sage -python src/solve7_main.py \
      --mode cuso \
      --a "$A" --b "$B" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph "$CUSO_GRAPH" \
      2>&1 | tee -a "$log"
    status=${PIPESTATUS[0]}
    echo "exit_status=$status" | tee -a "$log"
  ) &
  pids+=("$!")
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

grep -RniE "FOUND|p =|q =|m.hex|Traceback|Exception|cuso unavailable|cuso fail|not found|exit_status" "$LOG_DIR" || true
echo "logs written to $LOG_DIR"
