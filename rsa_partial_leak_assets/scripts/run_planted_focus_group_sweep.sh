#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_focus_group_sweep_$(date +%Y%m%d_%H%M%S)}"
PLANTED_BITS="${PLANTED_BITS:-64}"
PLANTED_SEED="${PLANTED_SEED:-20260705}"
LOW_FREE_WIDTH="${LOW_FREE_WIDTH:-1}"
CONSTRUCTION="${CONSTRUCTION:-y-only}"
X_SHIFT_LIMIT="${X_SHIFT_LIMIT:-1}"
X_ORIGIN="${X_ORIGIN:-zero}"
M_VALUES="${M_VALUES:-3-12}"
T_VALUES="${T_VALUES:-1-6}"
mkdir -p "$LOG_DIR"

python3 experiments/planted_focus_group_hm.py \
  --mode self-test \
  --prime-bits "$PLANTED_BITS" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  | tee "$LOG_DIR/self_test.txt"

python3 experiments/planted_focus_group_hm.py \
  --mode sweep \
  --prime-bits "$PLANTED_BITS" \
  --seed "$PLANTED_SEED" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m-values "$M_VALUES" \
  --t-values "$T_VALUES" \
  | tee "$LOG_DIR/sweep.txt"

python3 experiments/planted_focus_group_hm.py \
  --mode sweep \
  --prime-bits "$PLANTED_BITS" \
  --seed "$PLANTED_SEED" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m-values 12 \
  --t-values 3-4 \
  --json \
  | tee "$LOG_DIR/vanishing_candidates.json"

python3 experiments/planted_focus_group_hm.py \
  --mode prune \
  --prime-bits "$PLANTED_BITS" \
  --seed "$PLANTED_SEED" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m 12 \
  --t 3 \
  --json \
  | tee "$LOG_DIR/prune_m12_t3.json"

python3 experiments/planted_focus_group_hm.py \
  --mode prune \
  --prime-bits "$PLANTED_BITS" \
  --seed "$PLANTED_SEED" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m 12 \
  --t 4 \
  --json \
  | tee "$LOG_DIR/prune_m12_t4.json"

python3 experiments/planted_focus_group_hm.py \
  --mode drop-sweep \
  --prime-bits "$PLANTED_BITS" \
  --seed "$PLANTED_SEED" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m 12 \
  --t 3 \
  --drop-axes power,nscale \
  | tee "$LOG_DIR/drop_sweep_m12_t3.txt"

python3 experiments/planted_focus_group_hm.py \
  --mode drop-sweep \
  --prime-bits "$PLANTED_BITS" \
  --seed "$PLANTED_SEED" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m 12 \
  --t 4 \
  --drop-axes power,nscale \
  | tee "$LOG_DIR/drop_sweep_m12_t4.txt"

echo "construction=$CONSTRUCTION"
echo "x_origin=$X_ORIGIN"
echo "x_shift_limit=$X_SHIFT_LIMIT"
echo "logs written to $LOG_DIR"
