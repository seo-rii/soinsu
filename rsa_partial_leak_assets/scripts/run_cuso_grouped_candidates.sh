#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

A="${1:-0}"
B="${2:-16}"
LOG_DIR="${3:-logs/grouped_candidates_$(date +%Y%m%d_%H%M%S)}"

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"

CUSO_LOG="${CUSO_LOG:-INFO}"
CUSO_GRAPH="${CUSO_GRAPH:-off}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-60s}"
CUSO_UPPER_LOW_VARIABLES="${CUSO_UPPER_LOW_VARIABLES:-0}"
CUSO_UPPER_ALL_VARIABLES="${CUSO_UPPER_ALL_VARIABLES:-0}"
CUSO_UPPER_VARIABLES="${CUSO_UPPER_VARIABLES:-}"
STOP_ON_FOUND="${STOP_ON_FOUND:-1}"

mkdir -p "$LOG_DIR"

origin_args=()
if [[ "$CUSO_UPPER_LOW_VARIABLES" == "1" || "$CUSO_UPPER_LOW_VARIABLES" == "true" ]]; then
  origin_args+=(--cuso-upper-low-variables)
fi
if [[ "$CUSO_UPPER_ALL_VARIABLES" == "1" || "$CUSO_UPPER_ALL_VARIABLES" == "true" ]]; then
  origin_args+=(--cuso-upper-all-variables)
fi
for spec in $CUSO_UPPER_VARIABLES; do
  origin_args+=(--cuso-upper-variable "$spec")
done

if (( B <= A )); then
  echo "empty candidate range $A:$B" >&2
  exit 2
fi

for ((cid=A; cid<B; cid++)); do
  next=$((cid + 1))
  cid_tag="$(printf '%03d' "$cid")"
  log="$LOG_DIR/grouped_cid_${cid_tag}.log"
  echo "cid=$cid range=$cid:$next graph=$CUSO_GRAPH timeout=$CUSO_TIMEOUT upper_low=$CUSO_UPPER_LOW_VARIABLES upper_all=$CUSO_UPPER_ALL_VARIABLES upper_variables=$CUSO_UPPER_VARIABLES" | tee "$log"
  set +e
  timeout "$CUSO_TIMEOUT" sage -python src/solve7_main.py \
    --mode cuso \
    --a "$cid" --b "$next" \
    --cuso-log "$CUSO_LOG" \
    --cuso-graph "$CUSO_GRAPH" \
    "${origin_args[@]}" \
    2>&1 | tee -a "$log"
  status=${PIPESTATUS[0]}
  set -e
  echo "exit_status=$status" | tee -a "$log"
  if grep -q "FOUND" "$log"; then
    echo "factor found in cid=$cid" | tee -a "$LOG_DIR/FOUND.txt"
    if [[ "$STOP_ON_FOUND" == "1" || "$STOP_ON_FOUND" == "true" ]]; then
      break
    fi
  fi
done

python3 experiments/summarize_cuso_logs.py "$LOG_DIR" | tee "$LOG_DIR/summary.txt"
echo "candidate_range=$A:$B"
echo "logs written to $LOG_DIR"
