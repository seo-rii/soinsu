#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_lowfree_sweep_$(date +%Y%m%d_%H%M%S)}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-20s}"
PLANTED_BITS="${PLANTED_BITS:-64}"
PLANTED_SEED="${PLANTED_SEED:-20260705}"
LOW_FREE_WIDTHS="${LOW_FREE_WIDTHS:-0,1,2,4,8}"
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

run_step "bits${PLANTED_BITS}_list_widths" \
  python3 experiments/planted_lowfree_sweep.py \
    --mode list-widths \
    --prime-bits "$PLANTED_BITS" \
    --seed "$PLANTED_SEED" \
    --low-free-widths "$LOW_FREE_WIDTHS"

IFS=',' read -r -a widths <<< "$LOW_FREE_WIDTHS"
for width in "${widths[@]}"; do
  run_step "bits${PLANTED_BITS}_lowfree${width}_graph_off" \
    sage -python experiments/planted_lowfree_sweep.py \
      --mode cuso \
      --prime-bits "$PLANTED_BITS" \
      --seed "$PLANTED_SEED" \
      --low-free-widths "$width" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph off \
      --json
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "logs written to $LOG_DIR"
