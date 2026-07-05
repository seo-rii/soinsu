#!/usr/bin/env bash
set -euo pipefail
A="${1:-0}"
B="${2:-256}"
mkdir -p logs
sage -python src/solve7_main.py --mode cuso --a "$A" --b "$B" | tee "logs/cuso_${A}_${B}.log"
