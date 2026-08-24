#!/usr/bin/env python3
"""Summarize cuso smoke logs.

The cuso INFO logs are verbose and smoke runs often end by timeout.  This helper
extracts the last useful progress markers so option/shape runs can be compared
without reading each log manually.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


RE_EXIT = re.compile(r"exit_status=(\d+)")
RE_SHIFT = re.compile(r"Computed (\d+) shift polynomials")
RE_GRAPH = re.compile(r"Graph optimization found a subset of (\d+) shift relations \(out of (\d+)\)")
RE_LATTICE = re.compile(r"Built a dual lattice basis of rank (\d+) and dimension (\d+)")
RE_INTREL = re.compile(r"Found (\d+) integer relations")
RE_MULT = re.compile(r"Generated ideal of multiplicity (\d+) with ([0-9.]+)-bit modulus")
RE_ROOTS = re.compile(r"(?:cid \d+ roots|split roots) (\d+)")
RE_JSON_STATUS = re.compile(r'"status":\s*"([^"]+)"')
RE_JSON_ROOTS = re.compile(r'"roots":\s*(\d+)')
RE_QUEUE_INDEX = re.compile(r"queue_index=(\d+)")
RE_QUEUE_KEY = re.compile(r"(?:^|\s)key=([^\s]+)")
RE_QUEUE_SHAPE = re.compile(r"(?:^|\s)shape=([^\s]+)")


def summarize(path: Path) -> dict:
    out = {
        "log": str(path),
        "exit_status": None,
        "status": None,
        "max_multiplicity": None,
        "max_shift_polys": None,
        "graph_subset": None,
        "graph_total": None,
        "last_rank": None,
        "last_dimension": None,
        "last_integer_relations": None,
        "roots": None,
        "traceback": False,
        "queue_index": None,
        "queue_key": None,
        "queue_shape": None,
    }
    for line in path.read_text(errors="replace").splitlines():
        if "Traceback" in line or "Exception" in line:
            out["traceback"] = True
        if m := RE_EXIT.search(line):
            out["exit_status"] = int(m.group(1))
        if m := RE_JSON_STATUS.search(line):
            out["status"] = m.group(1)
        if m := RE_MULT.search(line):
            out["max_multiplicity"] = int(m.group(1))
        if m := RE_SHIFT.search(line):
            out["max_shift_polys"] = int(m.group(1))
        if m := RE_GRAPH.search(line):
            out["graph_subset"] = int(m.group(1))
            out["graph_total"] = int(m.group(2))
        if m := RE_LATTICE.search(line):
            out["last_rank"] = int(m.group(1))
            out["last_dimension"] = int(m.group(2))
        if m := RE_INTREL.search(line):
            out["last_integer_relations"] = int(m.group(1))
        if m := RE_ROOTS.search(line):
            out["roots"] = int(m.group(1))
        if m := RE_JSON_ROOTS.search(line):
            out["roots"] = int(m.group(1))
        if m := RE_QUEUE_INDEX.search(line):
            out["queue_index"] = int(m.group(1))
        if m := RE_QUEUE_KEY.search(line):
            out["queue_key"] = m.group(1)
        if m := RE_QUEUE_SHAPE.search(line):
            out["queue_shape"] = m.group(1)
        if "FOUND" in line:
            out["status"] = "factor"
    if out["status"] is None:
        if out["exit_status"] == 124:
            out["status"] = "timeout"
        elif out["traceback"]:
            out["status"] = "error"
        elif out["exit_status"] == 0:
            out["status"] = "ok"
        else:
            out["status"] = "incomplete"
    return out


def print_table(rows: list[dict]) -> None:
    headers = [
        "log",
        "queue",
        "key",
        "shape",
        "status",
        "exit",
        "mult",
        "shift",
        "graph",
        "rank",
        "intrel",
        "roots",
    ]
    data = []
    for row in rows:
        graph = ""
        if row["graph_subset"] is not None:
            graph = f"{row['graph_subset']}/{row['graph_total']}"
        rank = ""
        if row["last_rank"] is not None:
            rank = f"{row['last_rank']}/{row['last_dimension']}"
        data.append(
            [
                Path(row["log"]).name,
                "" if row["queue_index"] is None else str(row["queue_index"]),
                row["queue_key"] or "",
                row["queue_shape"] or "",
                row["status"] or "",
                "" if row["exit_status"] is None else str(row["exit_status"]),
                "" if row["max_multiplicity"] is None else str(row["max_multiplicity"]),
                "" if row["max_shift_polys"] is None else str(row["max_shift_polys"]),
                graph,
                rank,
                "" if row["last_integer_relations"] is None else str(row["last_integer_relations"]),
                "" if row["roots"] is None else str(row["roots"]),
            ]
        )
    widths = [len(h) for h in headers]
    for row in data:
        for i, value in enumerate(row):
            widths[i] = max(widths[i], len(value))
    print("  ".join(h.ljust(widths[i]) for i, h in enumerate(headers)))
    print("  ".join("-" * widths[i] for i in range(len(headers))))
    for row in data:
        print("  ".join(value.ljust(widths[i]) for i, value in enumerate(row)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    logs: list[Path] = []
    for raw in args.paths:
        path = Path(raw)
        if path.is_dir():
            logs.extend(sorted(path.rglob("*.log")))
        else:
            logs.append(path)
    rows = [summarize(path) for path in logs]
    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
    else:
        print_table(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
