#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/planted_cuso_scale_sweep_$(date +%Y%m%d_%H%M%S)}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-20s}"
PLANTED_SEED="${PLANTED_SEED:-20260705}"
PLANTED_BITS_LIST="${PLANTED_BITS_LIST:-64}"
PLANTED_SHAPES="${PLANTED_SHAPES:-low_tail_univariate S0_grouped_2var}"
PLANTED_OPTIONS="${PLANTED_OPTIONS:-graph_off graph_on graph_off_no_intermediate graph_off_allow_partial}"
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

for bits in $PLANTED_BITS_LIST; do
  run_step "bits${bits}_list_shapes" \
    python3 experiments/planted_cuso_smoke.py \
      --mode list-shapes \
      --prime-bits "$bits"

  for shape in $PLANTED_SHAPES; do
    for option in $PLANTED_OPTIONS; do
      args="$(option_args "$option")" || exit $?
      # shellcheck disable=SC2086
      run_step "bits${bits}_${shape}_${option}" \
        sage -python experiments/planted_cuso_smoke.py \
          --mode cuso \
          --shape "$shape" \
          --prime-bits "$bits" \
          --seed "$PLANTED_SEED" \
          --cuso-log "$CUSO_LOG" \
          $args \
          --json
    done
  done
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "logs written to $LOG_DIR"
