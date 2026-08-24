#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

LOG_DIR="${1:-logs/ct07_cuso/learned_leaf_cycle_$(date +%Y%m%d_%H%M%S)}"

CYCLES="${CYCLES:-1}"
SHAPES="${SHAPES:-F_150_4_265_84_362_50}"
PROJECTION_PHASES="${PROJECTION_PHASES:-random}"
PHASE_SEEDS="${PHASE_SEEDS:-89}"
ITERATIONS="${ITERATIONS:-1}"
PREFIX_BITS="${PREFIX_BITS:-600}"
CONF_BUDGET="${CONF_BUDGET:-100000}"
QUEUE_LIMIT="${QUEUE_LIMIT:-1}"
MAX_LEAF_BITS="${MAX_LEAF_BITS:-8}"
LEAF_SHARDS="${LEAF_SHARDS:-8}"
LEAF_TIMEOUT="${LEAF_TIMEOUT:-180s}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_LOG="${CUSO_LOG:-}"
CUSO_NO_INTERMEDIATE="${CUSO_NO_INTERMEDIATE:-0}"
CUSO_ALLOW_PARTIAL="${CUSO_ALLOW_PARTIAL:-0}"
INITIAL_LEARNED_JSONL="${INITIAL_LEARNED_JSONL:-}"
LEARNED_ALL_JSONL="${LEARNED_ALL_JSONL:-$LOG_DIR/learned_all.jsonl}"
SKIP_LEAF_SHARDS="${SKIP_LEAF_SHARDS:-0}"

mkdir -p "$LOG_DIR"

merge_inputs=()
for path in $INITIAL_LEARNED_JSONL; do
  merge_inputs+=("$path")
done
python3 experiments/merge_learned_clauses.py \
  "${merge_inputs[@]}" \
  --output-jsonl "$LEARNED_ALL_JSONL" \
  --ignore-missing \
  --json | tee "$LOG_DIR/learned_initial_summary.json"

if (( CYCLES < 1 )); then
  echo "cycles=$CYCLES; initialized $LEARNED_ALL_JSONL"
  exit 0
fi

for cycle in $(seq 1 "$CYCLES"); do
  cycle_dir="$LOG_DIR/cycle_${cycle}"
  collect_dir="$cycle_dir/collect"
  leaf_dir="$cycle_dir/leaf_shards"
  learned_cycle="$cycle_dir/learned.jsonl"
  learned_next="$cycle_dir/learned_all.next.jsonl"
  mkdir -p "$cycle_dir"

  echo "cycle=$cycle collect_dir=$collect_dir"
  PREFIX_BITS="$PREFIX_BITS" \
    SHAPES="$SHAPES" \
    PROJECTION_PHASES="$PROJECTION_PHASES" \
    PHASE_SEEDS="$PHASE_SEEDS" \
    ITERATIONS="$ITERATIONS" \
    CONF_BUDGET="$CONF_BUDGET" \
    LOAD_LEARNED_JSONL="$LEARNED_ALL_JSONL" \
    bash scripts/run_programmatic_low600_queue_collect.sh "$collect_dir"

  if [[ "$SKIP_LEAF_SHARDS" == "1" || "$SKIP_LEAF_SHARDS" == "true" ]]; then
    echo "cycle=$cycle skip_leaf_shards=1"
    continue
  fi

  echo "cycle=$cycle leaf_dir=$leaf_dir"
  QUEUE_LIMIT="$QUEUE_LIMIT" \
    MAX_LEAF_BITS="$MAX_LEAF_BITS" \
    LEAF_SHARDS="$LEAF_SHARDS" \
    LEAF_TIMEOUT="$LEAF_TIMEOUT" \
    CUSO_GRAPH="$CUSO_GRAPH" \
    CUSO_LOG="$CUSO_LOG" \
    CUSO_NO_INTERMEDIATE="$CUSO_NO_INTERMEDIATE" \
    CUSO_ALLOW_PARTIAL="$CUSO_ALLOW_PARTIAL" \
    bash scripts/run_programmatic_low600_queue_leaf_shards.sh \
      "$collect_dir/projections.jsonl" \
      "$leaf_dir"

  python3 experiments/export_leaf_learned.py \
    "$collect_dir/projections.jsonl" \
    "$leaf_dir" \
    --output-jsonl "$learned_cycle" \
    --json | tee "$cycle_dir/export_summary.json"

  python3 experiments/merge_learned_clauses.py \
    "$LEARNED_ALL_JSONL" \
    "$learned_cycle" \
    --output-jsonl "$learned_next" \
    --ignore-missing \
    --json | tee "$cycle_dir/learned_merge_summary.json"
  mv "$learned_next" "$LEARNED_ALL_JSONL"
done

echo "learned clauses: $LEARNED_ALL_JSONL"
echo "logs written to $LOG_DIR"
