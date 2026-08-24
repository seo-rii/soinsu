#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/cuso_smoke_matrix_$(date +%Y%m%d_%H%M%S)}"
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

run_step grouped_cid0_graph_on \
  sage -python src/solve7_main.py \
    --mode cuso \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on

run_step grouped_cid0_graph_off \
  sage -python src/solve7_main.py \
    --mode cuso \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off

run_step grouped_cid0_graph_on_no_intermediate \
  sage -python src/solve7_main.py \
    --mode cuso \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on \
    --cuso-no-intermediate

run_step split_edge0_graph_on \
  sage -python src/solve7_main.py \
    --mode cuso-split \
    --cuso-split-brute-small-edges \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on

run_step mixed_s2_low_exact_high_grouped_graph_off \
  sage -python experiments/mixed_shape_cuso.py \
    --mode cuso \
    --shape S2_low_exact_high_grouped \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step mixed_s3_low_grouped_high_exact_graph_off \
  sage -python experiments/mixed_shape_cuso.py \
    --mode cuso \
    --shape S3_low_grouped_high_exact \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step mixed_s4_low_exact_high_mixed_graph_off \
  sage -python experiments/mixed_shape_cuso.py \
    --mode cuso \
    --shape S4_low_exact_high_mixed \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step partial_shape_b_graph_on \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape B \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on \
    --json

run_step partial_shape_b_graph_off \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape B \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step partial_shape_c_graph_on \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape C \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on \
    --json

run_step partial_shape_c_graph_off \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape C \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step partial_shape_d_265_64_362_8_graph_off \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape D_265_64_362_8 \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step partial_shape_d_265_48_362_16_graph_off \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape D_265_48_362_16 \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    --json

run_step partial_shape_a_graph_on \
  sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape A \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on \
    --json

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "logs written to $LOG_DIR"
