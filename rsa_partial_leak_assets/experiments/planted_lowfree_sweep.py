#!/usr/bin/env python3
"""Planted cuso sweep for low-variable-width retention boundaries.

`planted_cuso_smoke.py` showed that a one-variable known-low/high-tail model
works at 64 bits while the grouped S0 broad model does not.  This experiment
opens only a suffix of the scaled low group as an additional variable and keeps
the rest fixed, so we can find where true-branch recovery starts failing.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path

ASSET_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ASSET_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from solve7_main import _patch_cuso_boundset_for_sage109, cuso_find_small_roots_compat  # noqa: E402
from planted_cuso_smoke import Segment, make_instance, scaled_segments  # noqa: E402


def synthetic_instance(prime_bits: int) -> tuple[int, int, int]:
    p = (1 << (prime_bits - 1)) | 0x123456789ABCDEF
    q = (1 << (prime_bits - 1)) | 0xFEDCBA987654321
    return p, q, p * q


@dataclass(frozen=True)
class Model:
    prime_bits: int
    low_free_width: int
    low_free_offset: int
    n: int
    p: int
    base: int
    variables: tuple[Segment, ...]
    fixed: tuple[tuple[Segment, int], ...]

    @property
    def variable_bits(self) -> int:
        return sum(seg.width for seg in self.variables)

    @property
    def fixed_bits(self) -> int:
        return sum(seg.width for seg, _ in self.fixed)


def variable_name(seg: Segment) -> str:
    return f"v{seg.start}_{seg.width}"


def build_model(
    prime_bits: int,
    p: int,
    n: int,
    low_free_width: int,
    low_free_offset: int | None = None,
) -> Model:
    segments = scaled_segments(prime_bits)
    low_a = segments["low_a"]
    low_b = segments["low_b"]
    high_a = segments["high_a"]
    edge_low = segments["edge_low"]
    low_group = Segment(low_a.start, low_b.stop - low_a.start, "low_group")
    high_tail = Segment(high_a.start, prime_bits - high_a.start, "high_tail")
    if not 0 <= low_free_width <= low_group.width:
        raise ValueError(f"low_free_width must be in 0..{low_group.width}")
    if low_free_offset is None:
        low_free_offset = low_group.width - low_free_width
    if low_free_width == 0:
        low_free_offset = 0
    if not 0 <= low_free_offset <= low_group.width - low_free_width:
        raise ValueError(
            "low_free_offset must allow the free window to fit inside low_group"
        )

    variables = [high_tail]
    fixed_segments = [edge_low, low_group]
    if low_free_width:
        low_free = Segment(low_group.start + low_free_offset, low_free_width, "low_free")
        variables.insert(0, low_free)
        fixed_segments = [edge_low]
        if low_free_offset:
            fixed_segments.append(Segment(low_group.start, low_free_offset, "low_fixed_left"))
        right_width = low_group.width - low_free_offset - low_free_width
        if right_width:
            fixed_segments.append(Segment(low_free.stop, right_width, "low_fixed_right"))

    base = p
    for seg in variables:
        base &= ~seg.mask
    fixed = []
    for seg in fixed_segments:
        if seg.width <= 0:
            continue
        value = (p >> seg.start) & ((1 << seg.width) - 1)
        base &= ~seg.mask
        base |= value << seg.start
        fixed.append((seg, value))
    return Model(
        prime_bits,
        low_free_width,
        low_free_offset,
        n,
        p,
        base,
        tuple(variables),
        tuple(fixed),
    )


def record(model: Model, status: str, roots_count: int | None = None, elapsed: float | None = None) -> dict:
    out = {
        "status": status,
        "prime_bits": model.prime_bits,
        "low_free_width": model.low_free_width,
        "low_free_offset": model.low_free_offset,
        "n_bits": model.n.bit_length(),
        "fixed_bits": model.fixed_bits,
        "variable_bits": model.variable_bits,
        "fixed": [
            {"name": seg.name, "start": seg.start, "width": seg.width, "value": value}
            for seg, value in model.fixed
        ],
        "variables": [
            {"name": variable_name(seg), "segment": seg.name, "start": seg.start, "width": seg.width}
            for seg in model.variables
        ],
    }
    if roots_count is not None:
        out["roots"] = roots_count
    if elapsed is not None:
        out["elapsed"] = round(elapsed, 3)
    return out


def print_record(data: dict, as_json: bool) -> None:
    if as_json:
        print(json.dumps(data, sort_keys=True), flush=True)
        return
    variables = ", ".join(f"{item['segment']}:{item['width']}" for item in data["variables"])
    fixed = ", ".join(f"{item['name']}:{item['width']}={item['value']:#x}" for item in data["fixed"])
    parts = [
        f"status={data['status']}",
        f"prime_bits={data['prime_bits']}",
        f"low_free_width={data['low_free_width']}",
        f"low_free_offset={data['low_free_offset']}",
        f"variable_bits={data['variable_bits']}",
        f"fixed_bits={data['fixed_bits']}",
        f"fixed=[{fixed}]",
        f"variables=[{variables}]",
    ]
    if "roots" in data:
        parts.append(f"roots={data['roots']}")
    if "elapsed" in data:
        parts.append(f"elapsed={data['elapsed']}")
    print(" ".join(parts), flush=True)


def run_cuso(model: Model, args) -> tuple[str, int, float]:
    from sage.all import PolynomialRing, ZZ
    import cuso

    _patch_cuso_boundset_for_sage109()
    names = [variable_name(seg) for seg in model.variables]
    R = PolynomialRing(ZZ, names)
    gens = R.gens()
    f = R(model.base)
    bounds = {}
    for gen, seg in zip(gens, model.variables):
        f += R(1 << seg.start) * gen
        bounds[gen] = (0, 1 << seg.width)

    graph = None if args.cuso_graph == "auto" else args.cuso_graph == "on"
    t0 = time.time()
    roots = cuso_find_small_roots_compat(
        cuso,
        [f],
        bounds,
        modulus="p",
        modulus_multiple=model.n,
        modulus_lower_bound=1 << (model.prime_bits - 1),
        modulus_upper_bound=1 << model.prime_bits,
        use_graph_optimization=graph,
        use_intermediate_sizes=not args.cuso_no_intermediate,
        allow_partial_solutions=args.cuso_allow_partial,
    )
    elapsed = time.time() - t0
    status = "missed_true_root"
    for root in roots:
        candidate = model.base
        if isinstance(root, dict) and "p" in root:
            candidate = int(root["p"])
        else:
            ok = True
            for gen, seg in zip(gens, model.variables):
                value = None
                for key in (gen, str(gen)):
                    try:
                        if key in root:
                            value = int(root[key])
                            break
                    except Exception:
                        pass
                if value is None:
                    ok = False
                    break
                candidate += value << seg.start
            if not ok:
                continue
        if 1 < candidate < model.n and model.n % candidate == 0:
            status = "factor"
            break
        status = "candidate"
    return status, len(roots), elapsed


def parse_widths(text: str, max_width: int) -> list[int]:
    if text == "all":
        return list(range(max_width + 1))
    out = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a_s, b_s = part.split("-", 1)
            out.extend(range(int(a_s), int(b_s) + 1))
        else:
            out.append(int(part))
    return sorted(set(width for width in out if 0 <= width <= max_width))


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["dry-run", "cuso", "self-test", "list-widths"], default="dry-run")
    ap.add_argument("--prime-bits", type=int, default=64)
    ap.add_argument("--seed", type=int, default=20260705)
    ap.add_argument("--low-free-widths", default="0,1,2,4,8")
    ap.add_argument(
        "--low-free-offset",
        type=int,
        help="offset from the scaled low_group start; default keeps the suffix",
    )
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--cuso-log", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--cuso-graph", choices=["auto", "on", "off"], default="off")
    ap.add_argument("--cuso-no-intermediate", action="store_true")
    ap.add_argument("--cuso-allow-partial", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.cuso_log:
        logging.basicConfig(
            level=getattr(logging, args.cuso_log),
            format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        )

    if args.mode == "cuso":
        p, _q, n = make_instance(args.prime_bits, args.seed)
    else:
        p, _q, n = synthetic_instance(args.prime_bits)
    low_group = build_model(args.prime_bits, p, n, 0).fixed[-1][0]
    widths = parse_widths(args.low_free_widths, low_group.width)

    if args.mode == "self-test":
        for width in widths:
            model = build_model(args.prime_bits, p, n, width, args.low_free_offset)
            rebuilt = model.base
            for seg in model.variables:
                rebuilt += ((p >> seg.start) & ((1 << seg.width) - 1)) << seg.start
            assert rebuilt == p, (width, hex(rebuilt), hex(p))
        print("self-test ok")
        return 0

    if args.mode == "list-widths":
        for width in widths:
            print_record(
                record(
                    build_model(args.prime_bits, p, n, width, args.low_free_offset),
                    "dry_run",
                ),
                args.json,
            )
        return 0

    exit_status = 0 if args.mode == "dry-run" else 1
    for width in widths:
        model = build_model(args.prime_bits, p, n, width, args.low_free_offset)
        if args.mode == "dry-run":
            print_record(record(model, "dry_run"), args.json)
            continue
        status, roots_count, elapsed = run_cuso(model, args)
        print_record(record(model, status, roots_count, elapsed), args.json)
        if status == "factor":
            exit_status = 0
        else:
            exit_status = 1
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
