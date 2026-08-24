#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_cuso_smoke}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-30s}"
PLANTED_BITS="${PLANTED_BITS:-128}"
PLANTED_SEED="${PLANTED_SEED:-20260705}"
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

run_step self_test \
  python3 experiments/planted_cuso_smoke.py \
    --mode self-test \
    --prime-bits "$PLANTED_BITS"

run_step list_shapes \
  python3 experiments/planted_cuso_smoke.py \
    --mode list-shapes \
    --prime-bits "$PLANTED_BITS"

for shape in low_tail_univariate S0_grouped_2var S2_low_exact_high_grouped partial_B; do
  run_step "${shape}_cuso" \
    run_with_timeout sage -python experiments/planted_cuso_smoke.py \
      --mode cuso \
      --shape "$shape" \
      --prime-bits "$PLANTED_BITS" \
      --seed "$PLANTED_SEED" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph "$CUSO_GRAPH" \
      --json
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "logs written to $LOG_DIR"
