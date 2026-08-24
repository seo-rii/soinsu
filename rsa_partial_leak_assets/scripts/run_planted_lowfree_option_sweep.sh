#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_lowfree_option_sweep_$(date +%Y%m%d_%H%M%S)}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-30s}"
PLANTED_BITS="${PLANTED_BITS:-64}"
PLANTED_SEED="${PLANTED_SEED:-20260705}"
LOW_FREE_WIDTHS="${LOW_FREE_WIDTHS:-1}"
LOWFREE_OPTIONS="${LOWFREE_OPTIONS:-graph_off graph_on graph_auto graph_off_no_intermediate graph_off_allow_partial}"
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

option_args() {
  local option="$1"
  case "$option" in
    graph_off)
      echo "--cuso-graph off"
      ;;
    graph_on)
      echo "--cuso-graph on"
      ;;
    graph_auto)
      echo "--cuso-graph auto"
      ;;
    graph_off_no_intermediate)
      echo "--cuso-graph off --cuso-no-intermediate"
      ;;
    graph_off_allow_partial)
      echo "--cuso-graph off --cuso-allow-partial"
      ;;
    *)
      echo "unknown option: $option" >&2
      return 2
      ;;
  esac
}

run_step "bits${PLANTED_BITS}_list_widths" \
  python3 experiments/planted_lowfree_sweep.py \
    --mode list-widths \
    --prime-bits "$PLANTED_BITS" \
    --seed "$PLANTED_SEED" \
    --low-free-widths "$LOW_FREE_WIDTHS"

for width in $LOW_FREE_WIDTHS; do
  for option in $LOWFREE_OPTIONS; do
    args="$(option_args "$option")" || exit $?
    # shellcheck disable=SC2086
    run_step "bits${PLANTED_BITS}_lowfree${width}_${option}" \
      sage -python experiments/planted_lowfree_sweep.py \
        --mode cuso \
        --prime-bits "$PLANTED_BITS" \
        --seed "$PLANTED_SEED" \
        --low-free-widths "$width" \
        --cuso-log "$CUSO_LOG" \
        $args \
        --json
  done
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "logs written to $LOG_DIR"
