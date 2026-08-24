#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.py-site:$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/programmatic_low600_sat_cas_smoke}"
PREFIX_BITS="${PREFIX_BITS:-600}"
ITERATIONS="${ITERATIONS:-2}"
CONF_BUDGET="${CONF_BUDGET:-100000}"
PROJECTION_PHASE="${PROJECTION_PHASE:-default}"
PHASE_SEED="${PHASE_SEED:-20260706}"
RUN_CUSO="${RUN_CUSO:-0}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-45s}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_LOG="${CUSO_LOG:-INFO}"
QUEUE_JSONL="${QUEUE_JSONL:-$LOG_DIR/projections.jsonl}"
RESULTS_JSONL="${RESULTS_JSONL:-$LOG_DIR/results.jsonl}"
REPLAY_LIMIT="${REPLAY_LIMIT:-4}"
CUSO_REPLAY_LIMIT="${CUSO_REPLAY_LIMIT:-1}"
mkdir -p "$LOG_DIR"
: > "$QUEUE_JSONL"
: > "$RESULTS_JSONL"

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
  python3 experiments/programmatic_low600_sat_cas.py --mode self-test

run_step prefix_shape_b_dry_run \
  python3 experiments/programmatic_low600_sat_cas.py \
    --mode prefix-loop \
    --prefix-bits "$PREFIX_BITS" \
    --shape B \
    --iterations "$ITERATIONS" \
    --conf-budget "$CONF_BUDGET" \
    --projection-phase "$PROJECTION_PHASE" \
    --phase-seed "$PHASE_SEED" \
    --oracle dry-run \
    --queue-jsonl "$QUEUE_JSONL" \
    --json

run_step prefix_shape_c_dry_run \
  python3 experiments/programmatic_low600_sat_cas.py \
    --mode prefix-loop \
    --prefix-bits "$PREFIX_BITS" \
    --shape C \
    --iterations "$ITERATIONS" \
    --conf-budget "$CONF_BUDGET" \
    --projection-phase "$PROJECTION_PHASE" \
    --phase-seed "$PHASE_SEED" \
    --oracle dry-run \
    --queue-jsonl "$QUEUE_JSONL" \
    --json

run_step replay_queue_dry_run \
  python3 experiments/programmatic_low600_sat_cas.py \
    --mode replay-queue \
    --queue-jsonl "$QUEUE_JSONL" \
    --queue-limit "$REPLAY_LIMIT" \
    --oracle dry-run \
    --results-jsonl "$RESULTS_JSONL" \
    --json

if [[ "$RUN_CUSO" == "1" || "$RUN_CUSO" == "true" ]]; then
  run_step replay_queue_cuso \
    timeout "$CUSO_TIMEOUT" sage -python experiments/programmatic_low600_sat_cas.py \
      --mode replay-queue \
      --queue-jsonl "$QUEUE_JSONL" \
      --queue-limit "$CUSO_REPLAY_LIMIT" \
      --oracle cuso \
      --cuso-graph "$CUSO_GRAPH" \
      --cuso-log "$CUSO_LOG" \
      --results-jsonl "$RESULTS_JSONL" \
      --json
fi

echo "logs written to $LOG_DIR"
echo "queue written to $QUEUE_JSONL"
echo "results written to $RESULTS_JSONL"
