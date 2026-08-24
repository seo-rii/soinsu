#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

QUEUE_JSONL="${1:?usage: scripts/run_programmatic_low600_queue_focus_template.sh QUEUE_JSONL [LOG_DIR]}"
LOG_DIR="${2:-logs/programmatic_queue_focus_template_$(date +%Y%m%d_%H%M%S)}"

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.py-site:$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

QUEUE_START="${QUEUE_START:-0}"
QUEUE_LIMIT="${QUEUE_LIMIT:-}"
FOCUS_TIMEOUT="${FOCUS_TIMEOUT:-20s}"
FOCUS_X_ORIGIN="${FOCUS_X_ORIGIN:-both}"
FOCUS_RECOVER="${FOCUS_RECOVER:-both}"
FOCUS_PRIME_COUNT="${FOCUS_PRIME_COUNT:-56}"
FOCUS_ROWS="${FOCUS_ROWS:-13}"
FOCUS_M="${FOCUS_M:-11}"
FOCUS_T="${FOCUS_T:-3}"
FOCUS_BRUTE_SMALL_VARS="${FOCUS_BRUTE_SMALL_VARS:-0}"
FOCUS_NO_STOP_ON_EMPTY_PRIME="${FOCUS_NO_STOP_ON_EMPTY_PRIME:-0}"
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

focus_args=()
if [[ "$FOCUS_NO_STOP_ON_EMPTY_PRIME" == "1" || "$FOCUS_NO_STOP_ON_EMPTY_PRIME" == "true" ]]; then
  focus_args+=(--focus-no-stop-on-empty-prime)
fi

for row in "${queue_rows[@]}"; do
  IFS=$'\t' read -r queue_index key shape fixed_id <<< "$row"
  safe_key="$(printf '%s' "$key" | tr -c 'A-Za-z0-9_.-' '_')"
  log="$LOG_DIR/queue_${queue_index}_${safe_key}.log"
  echo "queue_index=$queue_index key=$key shape=$shape fixed_id=$fixed_id timeout=$FOCUS_TIMEOUT brute_small_vars=$FOCUS_BRUTE_SMALL_VARS" | tee "$log"
  set +e
  timeout "$FOCUS_TIMEOUT" python3 experiments/programmatic_low600_sat_cas.py \
    --mode replay-queue \
    --queue-jsonl "$QUEUE_JSONL" \
    --queue-start "$queue_index" \
    --queue-limit 1 \
    --oracle focus-template \
    --focus-x-origin "$FOCUS_X_ORIGIN" \
    --focus-recover "$FOCUS_RECOVER" \
    --focus-prime-count "$FOCUS_PRIME_COUNT" \
    --focus-rows "$FOCUS_ROWS" \
    --focus-m "$FOCUS_M" \
    --focus-t "$FOCUS_T" \
    --focus-brute-small-vars "$FOCUS_BRUTE_SMALL_VARS" \
    --results-jsonl "$RESULTS_JSONL" \
    --json \
    "${focus_args[@]}" \
    2>&1 | tee -a "$log"
  status=${PIPESTATUS[0]}
  set -e
  echo "exit_status=$status" | tee -a "$log"
done

python3 experiments/summarize_projection_queue.py "$QUEUE_JSONL" --results-jsonl "$RESULTS_JSONL" | tee "$LOG_DIR/queue_summary.txt"
python3 - "$RESULTS_JSONL" <<'PY' | tee "$LOG_DIR/focus_summary.txt"
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

records = []
with Path(sys.argv[1]).open("r", encoding="utf-8") as src:
    for line in src:
        if not line.strip():
            continue
        record = json.loads(line)
        if record.get("record_type") == "oracle_result":
            records.append(record)

status = Counter(record.get("status", "?") for record in records)
empty = Counter(record.get("empty_prime", "-") for record in records)
shape = Counter(record.get("shape", "?") for record in records)
origin = Counter(record.get("x_origin", "?") for record in records)
recover_times = [float(record.get("recover_elapsed", 0.0)) for record in records]
candidates = sum(int(record.get("candidate_count", 0)) for record in records)
print(f"oracle_results={len(records)}")
print("status:", " ".join(f"{key}={value}" for key, value in sorted(status.items(), key=lambda item: str(item[0]))))
print("empty_prime:", " ".join(f"{key}={value}" for key, value in sorted(empty.items(), key=lambda item: str(item[0]))))
print("shape:", " ".join(f"{key}={value}" for key, value in sorted(shape.items(), key=lambda item: str(item[0]))))
print("origin:", " ".join(f"{key}={value}" for key, value in sorted(origin.items(), key=lambda item: str(item[0]))))
print(f"candidate_count_total={candidates}")
if recover_times:
    print(f"recover_avg={statistics.fmean(recover_times):.3f}")
    print(f"recover_max={max(recover_times):.3f}")
PY
echo "queue=$QUEUE_JSONL"
echo "results=$RESULTS_JSONL"
echo "logs written to $LOG_DIR"
