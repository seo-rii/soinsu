#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/smoke}"
CUSO_LOG="${CUSO_LOG:-INFO}"
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

run_with_timeout() {
  if [[ "$CUSO_TIMEOUT" == "0" || -z "$CUSO_TIMEOUT" ]]; then
    "$@"
  else
    timeout "$CUSO_TIMEOUT" "$@"
  fi
}

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

run_step analyze run_with_timeout sage -python src/solve7_main.py --mode analyze

run_step grouped_cid0_graph_on \
  run_with_timeout sage -python src/solve7_main.py \
    --mode cuso \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on \
    "${origin_args[@]}"

run_step grouped_cid0_graph_off \
  run_with_timeout sage -python src/solve7_main.py \
    --mode cuso \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph off \
    "${origin_args[@]}"

run_step split_edge0_graph_on \
  run_with_timeout sage -python src/solve7_main.py \
    --mode cuso-split \
    --cuso-split-brute-small-edges \
    --a 0 --b 1 \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph on \
    "${origin_args[@]}"

echo "logs written to $LOG_DIR"
