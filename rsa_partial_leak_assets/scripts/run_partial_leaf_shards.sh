#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

SHAPE="${1:?usage: scripts/run_partial_leaf_shards.sh SHAPE FIXED_ID [LOG_DIR]}"
FIXED_ID="${2:?usage: scripts/run_partial_leaf_shards.sh SHAPE FIXED_ID [LOG_DIR]}"
LOG_DIR="${3:-logs/partial_leaf_shards_$(date +%Y%m%d_%H%M%S)}"

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.py-site:$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

LEAF_TOTAL="${LEAF_TOTAL:-256}"
LEAF_SHARDS="${LEAF_SHARDS:-8}"
LEAF_TIMEOUT="${LEAF_TIMEOUT:-180s}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_LOG="${CUSO_LOG:-}"
CUSO_NO_INTERMEDIATE="${CUSO_NO_INTERMEDIATE:-0}"
CUSO_ALLOW_PARTIAL="${CUSO_ALLOW_PARTIAL:-0}"

mkdir -p "$LOG_DIR"

common_args=()
if [[ -n "$CUSO_LOG" ]]; then
  common_args+=(--cuso-log "$CUSO_LOG")
fi
if [[ "$CUSO_NO_INTERMEDIATE" == "1" || "$CUSO_NO_INTERMEDIATE" == "true" ]]; then
  common_args+=(--cuso-no-intermediate)
fi
if [[ "$CUSO_ALLOW_PARTIAL" == "1" || "$CUSO_ALLOW_PARTIAL" == "true" ]]; then
  common_args+=(--cuso-allow-partial)
fi

shard_size=$(((LEAF_TOTAL + LEAF_SHARDS - 1) / LEAF_SHARDS))
statuses=()

for shard in $(seq 0 $((LEAF_SHARDS - 1))); do
  start=$((shard * shard_size))
  stop=$((start + shard_size))
  if (( start >= LEAF_TOTAL )); then
    continue
  fi
  if (( stop > LEAF_TOTAL )); then
    stop=$LEAF_TOTAL
  fi
  log="$LOG_DIR/shard_${shard}_${start}_${stop}.log"
  results="$LOG_DIR/shard_${shard}_${start}_${stop}.jsonl"
  (
    echo "shape=$SHAPE fixed_id=$FIXED_ID shard=$shard range=$start:$stop timeout=$LEAF_TIMEOUT"
    timeout "$LEAF_TIMEOUT" sage -python experiments/low600_partial_cuso.py \
      --mode leaf-certify \
      --shape "$SHAPE" \
      --a "$FIXED_ID" \
      --leaf-start "$start" \
      --leaf-stop "$stop" \
      --leaf-max-completions "$shard_size" \
      --cuso-graph "$CUSO_GRAPH" \
      --results-jsonl "$results" \
      --json \
      "${common_args[@]}"
    status=$?
    echo "exit_status=$status"
    exit "$status"
  ) > "$log" 2>&1 &
  statuses+=("$!")
done

overall=0
for pid in "${statuses[@]}"; do
  if ! wait "$pid"; then
    code=$?
    if [[ "$code" != "1" ]]; then
      overall=$code
    fi
  fi
done

python3 - "$LOG_DIR" <<'PY' | tee "$LOG_DIR/summary.txt"
import json
import sys
from collections import Counter
from pathlib import Path

log_dir = Path(sys.argv[1])
records = []
for path in sorted(log_dir.glob("shard_*.jsonl")):
    with path.open("r", encoding="utf-8") as src:
        for line in src:
            line = line.strip()
            if line:
                record = json.loads(line)
                record["_file"] = path.name
                records.append(record)

leaf_records = [record for record in records if "leaf_completion_id" in record]
summary_records = [record for record in records if "leaf_completion_id" not in record]
leaf_status = Counter(record.get("status", "?") for record in leaf_records)
summary_status = Counter(record.get("status", "?") for record in summary_records)
covered = sorted({record["leaf_completion_id"] for record in leaf_records})
missing = []
if summary_records:
    total = max(record.get("leaf_completions", 0) for record in summary_records)
    missing = [idx for idx in range(total) if idx not in set(covered)]
print(f"leaf_records={len(leaf_records)} summaries={len(summary_records)}")
print("leaf_status:", " ".join(f"{key}={value}" for key, value in sorted(leaf_status.items())) or "(none)")
print("summary_status:", " ".join(f"{key}={value}" for key, value in sorted(summary_status.items())) or "(none)")
print(f"covered={len(covered)} missing={len(missing)}")
if missing:
    print("missing:", ",".join(map(str, missing[:64])))
PY

echo "logs written to $LOG_DIR"
exit "$overall"
