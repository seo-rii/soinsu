#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/cuso_option_sweep_$(date +%Y%m%d_%H%M%S)}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-60s}"
mkdir -p "$LOG_DIR"

run_step() {
  local name="$1"
  shift
  local log="$LOG_DIR/${name}.log"
  echo "== $name ==" | tee "$log"
  set +e
  timeout "$CUSO_TIMEOUT" "$@" 2>&1 | tee -a "$log"
  local status=${PIPESTATUS[0]}
  set -e
  echo "exit_status=$status" | tee -a "$log"
}

run_mixed() {
  local name="$1"
  local shape="$2"
  shift 2
  run_step "$name" \
    sage -python experiments/mixed_shape_cuso.py \
      --mode cuso \
      --shape "$shape" \
      --a 0 --b 1 \
      --cuso-log "$CUSO_LOG" \
      "$@" \
      --json
}

run_partial() {
  local name="$1"
  local shape="$2"
  shift 2
  run_step "$name" \
    sage -python experiments/low600_partial_cuso.py \
      --mode cuso \
      --shape "$shape" \
      --a 0 --b 1 \
      --cuso-log "$CUSO_LOG" \
      "$@" \
      --json
}

for shape in S0_grouped_2var S2_low_exact_high_grouped; do
  prefix="$(echo "$shape" | tr '[:upper:]' '[:lower:]')"
  run_mixed "${prefix}_graph_off" "$shape" --cuso-graph off
  run_mixed "${prefix}_graph_off_no_intermediate" "$shape" --cuso-graph off --cuso-no-intermediate
  run_mixed "${prefix}_graph_off_allow_partial" "$shape" --cuso-graph off --cuso-allow-partial
  run_mixed "${prefix}_graph_on" "$shape" --cuso-graph on
done

run_partial partial_b_graph_off B --cuso-graph off
run_partial partial_b_graph_off_no_intermediate B --cuso-graph off --cuso-no-intermediate
run_partial partial_b_graph_off_allow_partial B --cuso-graph off --cuso-allow-partial
run_partial partial_b_graph_on B --cuso-graph on

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "logs written to $LOG_DIR"
