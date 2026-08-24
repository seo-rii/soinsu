#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/partial_low600_smoke}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_GRAPH="${CUSO_GRAPH:-on}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-0}"
UPPER_LOW_VARIABLES="${UPPER_LOW_VARIABLES:-0}"
UPPER_ALL_VARIABLES="${UPPER_ALL_VARIABLES:-0}"
UPPER_VARIABLES="${UPPER_VARIABLES:-}"
mkdir -p "$LOG_DIR"

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

origin_args=()
if [[ "$UPPER_LOW_VARIABLES" == "1" ]]; then
  origin_args+=(--upper-low-variables)
fi
if [[ "$UPPER_ALL_VARIABLES" == "1" ]]; then
  origin_args+=(--upper-all-variables)
fi
if [[ -n "$UPPER_VARIABLES" ]]; then
  IFS=',' read -r -a upper_variable_items <<< "$UPPER_VARIABLES"
  for item in "${upper_variable_items[@]}"; do
    if [[ -n "$item" ]]; then
      origin_args+=(--upper-variable "$item")
    fi
  done
fi

run_step self_test \
  python3 experiments/low600_partial_cuso.py --mode self-test

run_step list_shapes \
  python3 experiments/low600_partial_cuso.py --mode list-shapes

run_step shape_b_dry_run \
  python3 experiments/low600_partial_cuso.py \
    --mode dry-run \
    --shape B \
    --a 0 --b 1 \
    "${origin_args[@]}"

run_step shape_b_cuso \
  run_with_timeout sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape B \
    --a 0 --b 1 \
    "${origin_args[@]}" \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph "$CUSO_GRAPH" \
    --json

run_step shape_c_cuso \
  run_with_timeout sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape C \
    --a 0 --b 1 \
    "${origin_args[@]}" \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph "$CUSO_GRAPH" \
    --json

run_step shape_d_265_64_362_8_cuso \
  run_with_timeout sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape D_265_64_362_8 \
    --a 0 --b 1 \
    "${origin_args[@]}" \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph "$CUSO_GRAPH" \
    --json

run_step shape_d_265_48_362_16_cuso \
  run_with_timeout sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape D_265_48_362_16 \
    --a 0 --b 1 \
    "${origin_args[@]}" \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph "$CUSO_GRAPH" \
    --json

run_step shape_a_cuso \
  run_with_timeout sage -python experiments/low600_partial_cuso.py \
    --mode cuso \
    --shape A \
    --a 0 --b 1 \
    "${origin_args[@]}" \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph "$CUSO_GRAPH" \
    --json

echo "upper_low_variables=$UPPER_LOW_VARIABLES"
echo "upper_all_variables=$UPPER_ALL_VARIABLES"
echo "upper_variables=$UPPER_VARIABLES"
echo "logs written to $LOG_DIR"
