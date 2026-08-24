#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_focus_group_scale_$(date +%Y%m%d_%H%M%S)}"
PLANTED_BITS_LIST="${PLANTED_BITS_LIST:-64 80 96}"
PLANTED_SEED="${PLANTED_SEED:-20260705}"
LOW_FREE_WIDTH="${LOW_FREE_WIDTH:-1}"
CONSTRUCTION="${CONSTRUCTION:-y-only}"
X_SHIFT_LIMIT="${X_SHIFT_LIMIT:-1}"
X_ORIGIN="${X_ORIGIN:-zero}"
M_VALUES="${M_VALUES:-12}"
T_VALUES="${T_VALUES:-3-4}"
PRUNE_M="${PRUNE_M:-12}"
PRUNE_T_VALUES="${PRUNE_T_VALUES:-3 4}"
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

for bits in $PLANTED_BITS_LIST; do
  run_step "bits${bits}_self_test" \
    python3 experiments/planted_focus_group_hm.py \
      --mode self-test \
      --prime-bits "$bits" \
      --low-free-width "$LOW_FREE_WIDTH" \
      --construction "$CONSTRUCTION" \
      --x-shift-limit "$X_SHIFT_LIMIT" \
      --x-origin "$X_ORIGIN"

  run_step "bits${bits}_sweep_m${M_VALUES}_t${T_VALUES}" \
    python3 experiments/planted_focus_group_hm.py \
      --mode sweep \
      --prime-bits "$bits" \
      --seed "$PLANTED_SEED" \
      --low-free-width "$LOW_FREE_WIDTH" \
      --construction "$CONSTRUCTION" \
      --x-shift-limit "$X_SHIFT_LIMIT" \
      --x-origin "$X_ORIGIN" \
      --m-values "$M_VALUES" \
      --t-values "$T_VALUES" \
      --json

  for t_value in $PRUNE_T_VALUES; do
    run_step "bits${bits}_prune_m${PRUNE_M}_t${t_value}" \
      python3 experiments/planted_focus_group_hm.py \
        --mode prune \
        --prime-bits "$bits" \
        --seed "$PLANTED_SEED" \
        --low-free-width "$LOW_FREE_WIDTH" \
        --construction "$CONSTRUCTION" \
        --x-shift-limit "$X_SHIFT_LIMIT" \
        --x-origin "$X_ORIGIN" \
        --m "$PRUNE_M" \
        --t "$t_value" \
        --json
  done
done

echo "construction=$CONSTRUCTION"
echo "x_origin=$X_ORIGIN"
echo "x_shift_limit=$X_SHIFT_LIMIT"
echo "logs written to $LOG_DIR"
