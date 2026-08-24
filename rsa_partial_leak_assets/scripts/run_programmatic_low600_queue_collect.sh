#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

export PYTHONPATH="$PWD/.py-site:$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LOG_DIR="${1:-logs/programmatic_queue_collect_$(date +%Y%m%d_%H%M%S)}"
QUEUE_JSONL="${QUEUE_JSONL:-$LOG_DIR/projections.jsonl}"
PREFIX_BITS="${PREFIX_BITS:-600}"
SHAPES="${SHAPES:-B C}"
PROJECTION_PHASES="${PROJECTION_PHASES:-random}"
PHASE_SEEDS="${PHASE_SEEDS:-1 3 7 11}"
ITERATIONS="${ITERATIONS:-2}"
CONF_BUDGET="${CONF_BUDGET:-100000}"
QUEUE_APPEND="${QUEUE_APPEND:-0}"
LOAD_LEARNED_JSONL="${LOAD_LEARNED_JSONL:-}"
mkdir -p "$LOG_DIR"
if [[ "$QUEUE_APPEND" != "1" && "$QUEUE_APPEND" != "true" ]]; then
  : > "$QUEUE_JSONL"
fi

learned_args=()
for path in $LOAD_LEARNED_JSONL; do
  learned_args+=(--load-learned-jsonl "$path")
done

run_collect() {
  local shape="$1"
  local phase="$2"
  local seed="$3"
  local log="$LOG_DIR/collect_${shape}_${phase}_${seed}.log"
  echo "shape=$shape phase=$phase seed=$seed prefix_bits=$PREFIX_BITS iterations=$ITERATIONS" | tee "$log"
  set +e
  python3 experiments/programmatic_low600_sat_cas.py \
    --mode prefix-loop \
    --prefix-bits "$PREFIX_BITS" \
    --shape "$shape" \
    --iterations "$ITERATIONS" \
    --conf-budget "$CONF_BUDGET" \
    --projection-phase "$phase" \
    --phase-seed "$seed" \
    --oracle dry-run \
    --queue-jsonl "$QUEUE_JSONL" \
    "${learned_args[@]}" \
    --json \
    2>&1 | tee -a "$log"
  local status=${PIPESTATUS[0]}
  set -e
  echo "exit_status=$status" | tee -a "$log"
}

for shape in $SHAPES; do
  for phase in $PROJECTION_PHASES; do
    for seed in $PHASE_SEEDS; do
      run_collect "$shape" "$phase" "$seed"
    done
  done
done

python3 experiments/summarize_projection_queue.py "$QUEUE_JSONL" | tee "$LOG_DIR/queue_summary.txt"
echo "queue written to $QUEUE_JSONL"
echo "logs written to $LOG_DIR"
