#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs
for i in $(seq 0 7); do
  A=$((i*32)); B=$(((i+1)*32))
  bash scripts/run_cuso_range.sh "$A" "$B" &
done
wait
