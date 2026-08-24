#!/usr/bin/env python3
"""Merge learned-clause JSONL files, de-duplicating projection keys."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path


def merge_learned_clauses(input_paths: list[str], output_jsonl: str, ignore_missing: bool) -> dict:
    by_key: dict[str, dict] = {}
    ordered_keys: list[str] = []
    total = 0
    duplicates = 0
    skipped_missing = 0
    skipped_non_learned = 0

    for input_text in input_paths:
        path = Path(input_text)
        if not path.exists():
            if ignore_missing:
                skipped_missing += 1
                continue
            raise FileNotFoundError(path)
        with path.open("r", encoding="utf-8") as src:
            for line in src:
                line = line.strip()
                if not line:
                    continue
                record = json.loads(line)
                if record.get("record_type") != "learned_projection_clause":
                    skipped_non_learned += 1
                    continue
                key = record.get("key")
                if key is None:
                    key = f"{record['shape']}:{record['fixed_id']}"
                    record["key"] = key
                total += 1
                if key in by_key:
                    duplicates += 1
                    continue
                by_key[key] = record
                ordered_keys.append(key)

    out_path = Path(output_jsonl)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        for key in ordered_keys:
            out.write(json.dumps(by_key[key], sort_keys=True) + "\n")

    return {
        "inputs": len(input_paths),
        "total": total,
        "written": len(ordered_keys),
        "duplicates": duplicates,
        "skipped_missing": skipped_missing,
        "skipped_non_learned": skipped_non_learned,
        "output": str(out_path),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="*", help="learned-clause JSONL files")
    ap.add_argument("--output-jsonl", required=True)
    ap.add_argument("--ignore-missing", action="store_true")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            a = tmp / "a.jsonl"
            b = tmp / "b.jsonl"
            out = tmp / "out.jsonl"
            a.write_text(
                json.dumps({"record_type": "learned_projection_clause", "key": "A:1", "v": 1}) + "\n"
                + json.dumps({"record_type": "learned_projection_clause", "key": "B:2", "v": 2}) + "\n",
                encoding="utf-8",
            )
            b.write_text(
                json.dumps({"record_type": "learned_projection_clause", "key": "A:1", "v": 9}) + "\n"
                + json.dumps({"record_type": "learned_projection_clause", "shape": "C", "fixed_id": 3, "v": 3})
                + "\n",
                encoding="utf-8",
            )
            summary = merge_learned_clauses([str(a), str(b), str(tmp / "missing.jsonl")], str(out), True)
            merged = [json.loads(line) for line in out.read_text(encoding="utf-8").splitlines()]
            assert summary["written"] == 3
            assert summary["duplicates"] == 1
            assert summary["skipped_missing"] == 1
            assert [record["key"] for record in merged[:2]] == ["A:1", "B:2"]
            assert merged[0]["v"] == 1
            assert merged[2]["key"] == "C:3"
        print(json.dumps({"self_test": "ok"}, sort_keys=True))
        return 0

    summary = merge_learned_clauses(args.inputs, args.output_jsonl, args.ignore_missing)
    if args.json:
        print(json.dumps(summary, sort_keys=True))
    else:
        print(
            " ".join(
                [
                    f"inputs={summary['inputs']}",
                    f"total={summary['total']}",
                    f"written={summary['written']}",
                    f"duplicates={summary['duplicates']}",
                    f"skipped_missing={summary['skipped_missing']}",
                    f"output={summary['output']}",
                ]
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
