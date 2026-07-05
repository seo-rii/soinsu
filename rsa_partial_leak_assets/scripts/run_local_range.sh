#!/usr/bin/env bash
set -euo pipefail
A="${1:-0}"
B="${2:-1}"
M="${3:-8}"
T="${4:-3}"
LEAD="${5:-y}"
mkdir -p logs
python3 src/solve7_main.py --mode local --m "$M" --t "$T" --lead "$LEAD" --a "$A" --b "$B" | tee "logs/local_${A}_${B}_m${M}_t${T}_${LEAD}.log"
