#!/usr/bin/env python3
"""Summarize programmatic low600 projection queues."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path


def load_jsonl(path: Path) -> list[dict]:
    records = []
    with path.open("r", encoding="utf-8") as src:
        for lineno, line in enumerate(src, 1):
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            record["_line"] = lineno
            records.append(record)
    return records


def summarize_queue(queue_records: list[dict], result_records: list[dict]) -> dict:
    projections = [record for record in queue_records if record.get("record_type") == "prefix_projection"]
    keys = [record.get("key", f"{record.get('shape')}:{record.get('fixed_id')}") for record in projections]
    key_counts = Counter(keys)
    by_shape = Counter(record.get("shape", "?") for record in projections)
    by_phase = Counter(record.get("projection_phase", "?") for record in projections)
    by_prefix = Counter(record.get("prefix_bits", "?") for record in projections)
    clause_lengths = Counter(record.get("blocking_clause_len", "?") for record in projections)

    result_status = Counter()
    result_by_key: dict[str, Counter] = defaultdict(Counter)
    for record in result_records:
        if record.get("record_type") != "oracle_result":
            continue
        status = record.get("status", "?")
        key = record.get("key", f"{record.get('shape')}:{record.get('fixed_id')}")
        result_status[status] += 1
        result_by_key[key][status] += 1

    duplicate_keys = {
        key: count
        for key, count in sorted(key_counts.items(), key=lambda item: (-item[1], item[0]))
        if count > 1
    }
    examples = []
    seen = set()
    for record in projections:
        key = record.get("key", f"{record.get('shape')}:{record.get('fixed_id')}")
        if key in seen:
            continue
        seen.add(key)
        examples.append(
            {
                "line": record["_line"],
                "key": key,
                "shape": record.get("shape"),
                "fixed_id": record.get("fixed_id"),
                "fixed_bits": record.get("fixed_bits"),
                "blocking_clause_len": record.get("blocking_clause_len"),
                "projection_phase": record.get("projection_phase"),
                "phase_seed": record.get("phase_seed"),
                "prefix_bits": record.get("prefix_bits"),
                "result_status": dict(result_by_key.get(key, {})),
            }
        )
        if len(examples) >= 10:
            break

    return {
        "records": len(projections),
        "unique_keys": len(key_counts),
        "duplicate_records": len(projections) - len(key_counts),
        "by_shape": dict(sorted(by_shape.items())),
        "by_phase": dict(sorted(by_phase.items())),
        "by_prefix": dict(sorted(by_prefix.items(), key=lambda item: str(item[0]))),
        "blocking_clause_len": dict(sorted(clause_lengths.items(), key=lambda item: str(item[0]))),
        "duplicate_keys": duplicate_keys,
        "result_status": dict(sorted(result_status.items())),
        "examples": examples,
    }


def print_table(summary: dict) -> None:
    print(f"records={summary['records']} unique_keys={summary['unique_keys']} duplicate_records={summary['duplicate_records']}")
    print("by_shape:", " ".join(f"{key}={value}" for key, value in summary["by_shape"].items()) or "(none)")
    print("by_phase:", " ".join(f"{key}={value}" for key, value in summary["by_phase"].items()) or "(none)")
    print("by_prefix:", " ".join(f"{key}={value}" for key, value in summary["by_prefix"].items()) or "(none)")
    print("clause_len:", " ".join(f"{key}={value}" for key, value in summary["blocking_clause_len"].items()) or "(none)")
    if summary["result_status"]:
        print("result_status:", " ".join(f"{key}={value}" for key, value in summary["result_status"].items()))
    if summary["duplicate_keys"]:
        print("duplicates:")
        for key, count in list(summary["duplicate_keys"].items())[:10]:
            print(f"  {key} x{count}")
    print("examples:")
    for record in summary["examples"]:
        result_status = record["result_status"] or {}
        result = ",".join(f"{key}:{value}" for key, value in sorted(result_status.items())) or "-"
        print(
            "  "
            f"line={record['line']} key={record['key']} shape={record['shape']} "
            f"fixed_bits={record['fixed_bits']} clause={record['blocking_clause_len']} "
            f"phase={record['projection_phase']} seed={record['phase_seed']} result={result}"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("queue_jsonl")
    ap.add_argument("--results-jsonl")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    queue_records = load_jsonl(Path(args.queue_jsonl))
    result_records = load_jsonl(Path(args.results_jsonl)) if args.results_jsonl else []
    summary = summarize_queue(queue_records, result_records)
    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print_table(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
