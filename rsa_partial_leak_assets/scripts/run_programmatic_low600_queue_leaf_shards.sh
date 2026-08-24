#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

QUEUE_JSONL="${1:?usage: scripts/run_programmatic_low600_queue_leaf_shards.sh QUEUE_JSONL [LOG_DIR]}"
LOG_DIR="${2:-logs/programmatic_queue_leaf_shards_$(date +%Y%m%d_%H%M%S)}"

QUEUE_START="${QUEUE_START:-0}"
QUEUE_LIMIT="${QUEUE_LIMIT:-}"
MAX_LEAF_BITS="${MAX_LEAF_BITS:-8}"
LEAF_SHARDS="${LEAF_SHARDS:-8}"
LEAF_TIMEOUT="${LEAF_TIMEOUT:-180s}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"

mkdir -p "$LOG_DIR"

mapfile -t queue_rows < <(
  python3 - "$QUEUE_JSONL" "$QUEUE_START" "$QUEUE_LIMIT" "$MAX_LEAF_BITS" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
start = int(sys.argv[2])
limit = None if sys.argv[3] == "" else int(sys.argv[3])
max_leaf_bits = int(sys.argv[4])
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
        fixed_bits = int(record.get("fixed_bits", record.get("blocking_clause_len", 0)))
        leaf_bits = 146 - fixed_bits
        if leaf_bits < 0:
            print(f"skip\t{index}\tnegative_leaf_bits\t{record.get('key', '')}", file=sys.stderr)
            continue
        if leaf_bits > max_leaf_bits:
            print(f"skip\t{index}\tleaf_bits_{leaf_bits}_gt_{max_leaf_bits}\t{record.get('key', '')}", file=sys.stderr)
            continue
        key = record.get("key", f"{record['shape']}:{record['fixed_id']}")
        leaf_total = 1 << leaf_bits
        print(f"{index}\t{key}\t{record['shape']}\t{record['fixed_id']}\t{fixed_bits}\t{leaf_total}")
        emitted += 1
PY
)

if [[ "${#queue_rows[@]}" == "0" ]]; then
  echo "no queue records selected from $QUEUE_JSONL" >&2
  exit 2
fi

overall=0
for row in "${queue_rows[@]}"; do
  IFS=$'\t' read -r queue_index key shape fixed_id fixed_bits leaf_total <<< "$row"
  safe_key="$(printf '%s' "$key" | tr -c 'A-Za-z0-9_.-' '_')"
  record_dir="$LOG_DIR/queue_${queue_index}_${safe_key}"
  echo "queue_index=$queue_index key=$key shape=$shape fixed_id=$fixed_id fixed_bits=$fixed_bits leaf_total=$leaf_total"
  LEAF_TOTAL="$leaf_total" \
    LEAF_SHARDS="$LEAF_SHARDS" \
    LEAF_TIMEOUT="$LEAF_TIMEOUT" \
    CUSO_GRAPH="$CUSO_GRAPH" \
    scripts/run_partial_leaf_shards.sh "$shape" "$fixed_id" "$record_dir"
  status=$?
  if [[ "$status" != "0" ]]; then
    overall=$status
  fi
done

python3 - "$LOG_DIR" <<'PY' | tee "$LOG_DIR/summary.txt"
import json
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
summaries = []
leaf_records = 0
for path in sorted(root.glob("queue_*/summary.txt")):
    data = {"path": str(path.parent)}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("leaf_records="):
            for part in line.split():
                key, value = part.split("=", 1)
                data[key] = int(value)
        elif line.startswith("leaf_status:"):
            data["leaf_status_line"] = line.split(":", 1)[1].strip()
        elif line.startswith("summary_status:"):
            data["summary_status_line"] = line.split(":", 1)[1].strip()
        elif line.startswith("covered="):
            for part in line.split():
                key, value = part.split("=", 1)
                data[key] = int(value)
    summaries.append(data)
    leaf_records += data.get("leaf_records", 0)

status_counter = Counter()
for jsonl in sorted(root.glob("queue_*/*.jsonl")):
    with jsonl.open("r", encoding="utf-8") as src:
        for line in src:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if "leaf_completion_id" not in record:
                status_counter[record.get("status", "?")] += 1

print(f"queue_records={len(summaries)} leaf_records={leaf_records}")
print("summary_status:", " ".join(f"{key}={value}" for key, value in sorted(status_counter.items())) or "(none)")
for item in summaries:
    print(
        f"{item['path']} leaf_records={item.get('leaf_records', 0)} "
        f"covered={item.get('covered', 0)} missing={item.get('missing', 0)} "
        f"leaf_status={item.get('leaf_status_line', '')} "
        f"summary_status={item.get('summary_status_line', '')}"
    )
PY

echo "queue=$QUEUE_JSONL"
echo "logs written to $LOG_DIR"
exit "$overall"
