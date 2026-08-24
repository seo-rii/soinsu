#!/usr/bin/env python3
"""Export verified leaf-shard parent projections as learned-clause records."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def load_queue(path: Path) -> dict[str, dict]:
    out = {}
    with path.open("r", encoding="utf-8") as src:
        for line in src:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get("record_type") != "prefix_projection":
                continue
            key = record.get("key", f"{record['shape']}:{record['fixed_id']}")
            out[key] = record
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("queue_jsonl")
    ap.add_argument("leaf_log_root")
    ap.add_argument("--output-jsonl", required=True)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    queue = load_queue(Path(args.queue_jsonl))
    learned = []
    skipped = []
    for record_dir in sorted(Path(args.leaf_log_root).glob("queue_*")):
        records = []
        for jsonl in sorted(record_dir.glob("shard_*.jsonl")):
            with jsonl.open("r", encoding="utf-8") as src:
                for line in src:
                    line = line.strip()
                    if line:
                        records.append(json.loads(line))
        leaf_records = [record for record in records if "leaf_completion_id" in record]
        summaries = [record for record in records if "leaf_completion_id" not in record]
        if not summaries:
            skipped.append({"path": str(record_dir), "reason": "no_summary"})
            continue
        parent_shape = summaries[0]["shape"]
        parent_fixed_id = summaries[0]["fixed_id"]
        key = f"{parent_shape}:{parent_fixed_id}"
        projection = queue.get(key)
        if projection is None:
            skipped.append({"path": str(record_dir), "key": key, "reason": "missing_queue_record"})
            continue
        total = max(record.get("leaf_completions", 0) for record in summaries)
        covered = {record["leaf_completion_id"] for record in leaf_records}
        missing = [idx for idx in range(total) if idx not in covered]
        leaf_status = Counter(record.get("status", "?") for record in leaf_records)
        summary_status = Counter(record.get("status", "?") for record in summaries)
        if missing:
            skipped.append({"path": str(record_dir), "key": key, "reason": "missing_leaves", "missing": len(missing)})
            continue
        if set(leaf_status) != {"soft_no_root"}:
            skipped.append({"path": str(record_dir), "key": key, "reason": "non_soft_leaf", "leaf_status": dict(leaf_status)})
            continue
        if set(summary_status) != {"leaf_range_exhausted_soft_no_root"}:
            skipped.append(
                {
                    "path": str(record_dir),
                    "key": key,
                    "reason": "non_exhausted_summary",
                    "summary_status": dict(summary_status),
                }
            )
            continue
        learned.append(
            {
                "record_type": "learned_projection_clause",
                "version": 1,
                "key": key,
                "shape": parent_shape,
                "fixed_id": parent_fixed_id,
                "fixed": projection["fixed"],
                "fixed_bits": projection["fixed_bits"],
                "blocking_clause_len": projection["blocking_clause_len"],
                "source_oracle": "leaf-shards",
                "source_status": "all_leaves_soft_no_root",
                "leaf_records": len(leaf_records),
                "leaf_completions": total,
                "summary_records": len(summaries),
                "source_path": str(record_dir),
            }
        )

    out_path = Path(args.output_jsonl)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        for record in learned:
            out.write(json.dumps(record, sort_keys=True) + "\n")

    summary = {
        "learned": len(learned),
        "skipped": len(skipped),
        "output": str(out_path),
        "skipped_records": skipped[:10],
    }
    if args.json:
        print(json.dumps(summary, sort_keys=True))
    else:
        print(f"learned={len(learned)} skipped={len(skipped)} output={out_path}")
        for record in learned[:10]:
            print(f"  {record['key']} clause={record['blocking_clause_len']} leaves={record['leaf_records']}")
        for record in skipped[:10]:
            print(f"  skipped {record}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
