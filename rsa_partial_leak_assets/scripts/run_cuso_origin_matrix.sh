#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/cuso_origin_matrix_$(date +%Y%m%d_%H%M%S)}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-30s}"
EDGE_A="${EDGE_A:-0}"
EDGE_B="${EDGE_B:-1}"
MIXED_SHAPES="${MIXED_SHAPES:-S0_grouped_2var S2_low_exact_high_grouped}"
PARTIAL_SHAPES="${PARTIAL_SHAPES:-B C}"
GRAPHS="${GRAPHS:-off on}"
ORIGINS="${ORIGINS:-lower upper-low}"
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

origin_args_for() {
  local origin="$1"
  case "$origin" in
    lower)
      ;;
    upper-low)
      printf '%s\n' --upper-low-variables
      ;;
    upper-all)
      printf '%s\n' --upper-all-variables
      ;;
    *)
      echo "unknown origin: $origin" >&2
      return 2
      ;;
  esac
}

run_mixed_dry() {
  local shape="$1"
  local origin="$2"
  local origin_tag="${origin//-/_}"
  local shape_tag
  shape_tag="$(echo "$shape" | tr '[:upper:]' '[:lower:]')"
  mapfile -t origin_args < <(origin_args_for "$origin")
  run_step "mixed_${shape_tag}_${origin_tag}_dry_run" \
    python3 experiments/mixed_shape_cuso.py \
      --mode dry-run \
      --shape "$shape" \
      --a "$EDGE_A" --b "$EDGE_B" \
      "${origin_args[@]}" \
      --json
}

run_mixed_cuso() {
  local shape="$1"
  local origin="$2"
  local graph="$3"
  local origin_tag="${origin//-/_}"
  local shape_tag
  shape_tag="$(echo "$shape" | tr '[:upper:]' '[:lower:]')"
  mapfile -t origin_args < <(origin_args_for "$origin")
  run_step "mixed_${shape_tag}_${origin_tag}_graph_${graph}" \
    sage -python experiments/mixed_shape_cuso.py \
      --mode cuso \
      --shape "$shape" \
      --a "$EDGE_A" --b "$EDGE_B" \
      "${origin_args[@]}" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph "$graph" \
      --json
}

run_partial_dry() {
  local shape="$1"
  local origin="$2"
  local origin_tag="${origin//-/_}"
  local shape_tag
  shape_tag="$(echo "$shape" | tr '[:upper:]' '[:lower:]')"
  mapfile -t origin_args < <(origin_args_for "$origin")
  run_step "partial_${shape_tag}_${origin_tag}_dry_run" \
    python3 experiments/low600_partial_cuso.py \
      --mode dry-run \
      --shape "$shape" \
      --a "$EDGE_A" --b "$EDGE_B" \
      "${origin_args[@]}" \
      --json
}

run_partial_cuso() {
  local shape="$1"
  local origin="$2"
  local graph="$3"
  local origin_tag="${origin//-/_}"
  local shape_tag
  shape_tag="$(echo "$shape" | tr '[:upper:]' '[:lower:]')"
  mapfile -t origin_args < <(origin_args_for "$origin")
  run_step "partial_${shape_tag}_${origin_tag}_graph_${graph}" \
    sage -python experiments/low600_partial_cuso.py \
      --mode cuso \
      --shape "$shape" \
      --a "$EDGE_A" --b "$EDGE_B" \
      "${origin_args[@]}" \
      --cuso-log "$CUSO_LOG" \
      --cuso-graph "$graph" \
      --json
}

for origin in $ORIGINS; do
  for shape in $MIXED_SHAPES; do
    run_mixed_dry "$shape" "$origin"
    for graph in $GRAPHS; do
      run_mixed_cuso "$shape" "$origin" "$graph"
    done
  done

  for shape in $PARTIAL_SHAPES; do
    run_partial_dry "$shape" "$origin"
    for graph in $GRAPHS; do
      run_partial_cuso "$shape" "$origin" "$graph"
    done
  done
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "mixed_shapes=$MIXED_SHAPES"
echo "partial_shapes=$PARTIAL_SHAPES"
echo "origins=$ORIGINS"
echo "graphs=$GRAPHS"
echo "edge_range=$EDGE_A:$EDGE_B"
echo "logs written to $LOG_DIR"
