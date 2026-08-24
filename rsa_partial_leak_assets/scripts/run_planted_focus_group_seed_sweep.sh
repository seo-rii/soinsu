#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_focus_group_seed_$(date +%Y%m%d_%H%M%S)}"
PLANTED_BITS="${PLANTED_BITS:-96}"
PLANTED_SEEDS="${PLANTED_SEEDS:-20260705-20260714}"
LOW_FREE_WIDTH="${LOW_FREE_WIDTH:-1}"
CONSTRUCTION="${CONSTRUCTION:-x-limited}"
X_SHIFT_LIMIT="${X_SHIFT_LIMIT:-1}"
X_ORIGIN="${X_ORIGIN:-zero}"
M_VALUES="${M_VALUES:-11}"
T_VALUES="${T_VALUES:-3}"
JSON="${JSON:-0}"
TEMPLATE_FAMILIES="${TEMPLATE_FAMILIES:-}"
mkdir -p "$LOG_DIR"

json_args=()
log_name="seed_sweep.txt"
if [[ "$JSON" == "1" ]]; then
  json_args=(--json)
  log_name="seed_sweep.json"
fi
template_args=()
if [[ -n "$TEMPLATE_FAMILIES" ]]; then
  template_args=(--template-families "$TEMPLATE_FAMILIES")
fi

python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep \
  --prime-bits "$PLANTED_BITS" \
  --seeds "$PLANTED_SEEDS" \
  --low-free-width "$LOW_FREE_WIDTH" \
  --construction "$CONSTRUCTION" \
  --x-shift-limit "$X_SHIFT_LIMIT" \
  --x-origin "$X_ORIGIN" \
  --m-values "$M_VALUES" \
  --t-values "$T_VALUES" \
  "${template_args[@]}" \
  "${json_args[@]}" \
  | tee "$LOG_DIR/$log_name"

echo "construction=$CONSTRUCTION"
echo "x_origin=$X_ORIGIN"
echo "x_shift_limit=$X_SHIFT_LIMIT"
echo "template_families_count=$(tr ',' '\n' <<<"$TEMPLATE_FAMILIES" | sed '/^$/d' | wc -l)"
echo "json=$JSON"
echo "logs written to $LOG_DIR"
