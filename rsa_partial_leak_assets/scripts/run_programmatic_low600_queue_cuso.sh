#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

QUEUE_JSONL="${1:?usage: scripts/run_programmatic_low600_queue_cuso.sh QUEUE_JSONL [LOG_DIR]}"
LOG_DIR="${2:-logs/programmatic_queue_cuso_$(date +%Y%m%d_%H%M%S)}"

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.py-site:$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

QUEUE_START="${QUEUE_START:-0}"
QUEUE_LIMIT="${QUEUE_LIMIT:-}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-60s}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_NO_INTERMEDIATE="${CUSO_NO_INTERMEDIATE:-0}"
CUSO_ALLOW_PARTIAL="${CUSO_ALLOW_PARTIAL:-0}"
UPPER_LOW_VARIABLES="${UPPER_LOW_VARIABLES:-0}"
UPPER_ALL_VARIABLES="${UPPER_ALL_VARIABLES:-0}"
UPPER_VARIABLES="${UPPER_VARIABLES:-}"
RESULTS_JSONL="${RESULTS_JSONL:-$LOG_DIR/results.jsonl}"
mkdir -p "$LOG_DIR"
: > "$RESULTS_JSONL"

mapfile -t queue_rows < <(
  python3 - "$QUEUE_JSONL" "$QUEUE_START" "$QUEUE_LIMIT" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
start = int(sys.argv[2])
limit = None if sys.argv[3] == "" else int(sys.argv[3])
emitted = 0
with path.open("r", encoding="utf-8") as src:
    for index, line in enumerate(src):
        if index < start:
            continue
        if limit is not None and emitted >= limit:
            break
        line = line.strip()
        if not line:
            continue
        record = json.loads(line)
        if record.get("record_type") != "prefix_projection":
            continue
        key = record.get("key", f"{record['shape']}:{record['fixed_id']}")
        print(f"{index}\t{key}\t{record['shape']}\t{record['fixed_id']}")
        emitted += 1
PY
)

if [[ "${#queue_rows[@]}" == "0" ]]; then
  echo "no queue records selected from $QUEUE_JSONL" >&2
  exit 2
fi

common_args=()
if [[ "$CUSO_NO_INTERMEDIATE" == "1" || "$CUSO_NO_INTERMEDIATE" == "true" ]]; then
  common_args+=(--cuso-no-intermediate)
fi
if [[ "$CUSO_ALLOW_PARTIAL" == "1" || "$CUSO_ALLOW_PARTIAL" == "true" ]]; then
  common_args+=(--cuso-allow-partial)
fi
if [[ "$UPPER_LOW_VARIABLES" == "1" || "$UPPER_LOW_VARIABLES" == "true" ]]; then
  common_args+=(--upper-low-variables)
fi
if [[ "$UPPER_ALL_VARIABLES" == "1" || "$UPPER_ALL_VARIABLES" == "true" ]]; then
  common_args+=(--upper-all-variables)
fi
for spec in $UPPER_VARIABLES; do
  common_args+=(--upper-variable "$spec")
done

for row in "${queue_rows[@]}"; do
  IFS=$'\t' read -r queue_index key shape fixed_id <<< "$row"
  safe_key="$(printf '%s' "$key" | tr -c 'A-Za-z0-9_.-' '_')"
  log="$LOG_DIR/queue_${queue_index}_${safe_key}.log"
  echo "queue_index=$queue_index key=$key shape=$shape fixed_id=$fixed_id graph=$CUSO_GRAPH timeout=$CUSO_TIMEOUT" | tee "$log"
  set +e
  timeout "$CUSO_TIMEOUT" sage -python experiments/programmatic_low600_sat_cas.py \
    --mode replay-queue \
    --queue-jsonl "$QUEUE_JSONL" \
    --queue-start "$queue_index" \
    --queue-limit 1 \
    --oracle cuso \
    --cuso-graph "$CUSO_GRAPH" \
    --cuso-log "$CUSO_LOG" \
    --results-jsonl "$RESULTS_JSONL" \
    --json \
    "${common_args[@]}" \
    2>&1 | tee -a "$log"
  status=${PIPESTATUS[0]}
  set -e
  echo "exit_status=$status" | tee -a "$log"
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "queue=$QUEUE_JSONL"
echo "results=$RESULTS_JSONL"
echo "logs written to $LOG_DIR"
