#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

LOG_DIR="${1:-logs/smoke}"
CUSO_LOG="${CUSO_LOG:-INFO}"
mkdir -p "$LOG_DIR"

run_step() {
  local name="$1"
  shift
  local log="$LOG_DIR/${name}.log"
  echo "== $name ==" | tee "$log"
  set +e
  "$@" 2>&1 | tee -a "$log"
  local status=${PIPESTATUS[0]}
  set -e
  echo "exit_status=$status" | tee -a "$log"
}

run_step analyze sage -python src/solve7_main.py --mode analyze

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

run_step split_edge0_graph_on \
  sage -python src/solve7_main.py \
    --mode cuso-split \
    --cuso-split-brute-small-edges \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on

echo "logs written to $LOG_DIR"
