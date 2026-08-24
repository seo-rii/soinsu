#!/usr/bin/env bash
set -euo pipefail
A="${1:-0}"
B="${2:-256}"
export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"
CUSO_TIMEOUT="${CUSO_TIMEOUT:-0}"
CUSO_UPPER_LOW_VARIABLES="${CUSO_UPPER_LOW_VARIABLES:-0}"
CUSO_UPPER_ALL_VARIABLES="${CUSO_UPPER_ALL_VARIABLES:-0}"
CUSO_UPPER_VARIABLES="${CUSO_UPPER_VARIABLES:-}"
mkdir -p logs
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
cmd=(sage -python src/solve7_main.py --mode cuso --a "$A" --b "$B" "${origin_args[@]}")
if [[ "$CUSO_TIMEOUT" == "0" || -z "$CUSO_TIMEOUT" ]]; then
  "${cmd[@]}" | tee "logs/cuso_${A}_${B}.log"
else
  timeout "$CUSO_TIMEOUT" "${cmd[@]}" | tee "logs/cuso_${A}_${B}.log"
fi
