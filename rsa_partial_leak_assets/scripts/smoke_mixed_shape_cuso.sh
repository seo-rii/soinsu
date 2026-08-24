#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/mixed_shape_smoke}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-60s}"
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
  python3 experiments/mixed_shape_cuso.py --mode self-test

run_step list_shapes \
  python3 experiments/mixed_shape_cuso.py --mode list-shapes "${origin_args[@]}"

for shape in \
  S0_grouped_2var \
  S1_exact_5var \
  S2_low_exact_high_grouped \
  S3_low_grouped_high_exact \
  S4_low_exact_high_mixed
do
  run_step "${shape}_dry_run" \
    python3 experiments/mixed_shape_cuso.py \
      --mode dry-run \
      --shape "$shape" \
      --a 0 --b 1 \
      "${origin_args[@]}"

  run_step "${shape}_cuso" \
    run_with_timeout sage -python experiments/mixed_shape_cuso.py \
      --mode cuso \
      --shape "$shape" \
      --a 0 --b 1 \
      "${origin_args[@]}" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph "$CUSO_GRAPH" \
      --json
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "upper_low_variables=$UPPER_LOW_VARIABLES"
echo "upper_all_variables=$UPPER_ALL_VARIABLES"
echo "upper_variables=$UPPER_VARIABLES"
echo "logs written to $LOG_DIR"
