#!/usr/bin/env python3
"""Downscaled planted cuso smoke for true-branch retention checks.

This script creates a small RSA instance with the same relative unknown-block
layout as the challenge.  It then runs cuso on the true fixed branch for S0,
S2, or partial-B style models.  A successful planted run is not evidence that
the full challenge is easy, but it is a useful soundness gate before turning
heuristic no-root observations into SAT clauses.
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path

ASSET_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ASSET_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from solve7_main import (  # noqa: E402
    _patch_cuso_boundset_for_sage109,
    cuso_find_small_roots_compat,
)


@dataclass(frozen=True, order=True)
class Segment:
    start: int
    width: int
    name: str

    @property
    def stop(self) -> int:
        return self.start + self.width

    @property
    def mask(self) -> int:
        return ((1 << self.width) - 1) << self.start

    def label(self) -> str:
        return f"{self.name}:{self.start}:{self.width}"


ORIGINAL_SEGMENTS = (
    ("edge_low", 150, 4),
    ("low_a", 265, 84),
    ("low_b", 362, 58),
    ("high_a", 600, 69),
    ("high_b", 682, 87),
    ("high_c", 784, 46),
    ("edge_high", 920, 4),
)
SHAPES = ("S0_grouped_2var", "S2_low_exact_high_grouped", "partial_B", "low_tail_univariate")


def scaled_segments(prime_bits: int) -> dict[str, Segment]:
    if prime_bits < 64:
        raise ValueError("prime_bits must be at least 64 for this scaled geometry")
    out: dict[str, Segment] = {}
    cursor = 0
    for name, original_start, original_width in ORIGINAL_SEGMENTS:
        start = (original_start * prime_bits) // 1024
        width = max(1, (original_width * prime_bits + 512) // 1024)
        if start <= cursor:
            start = cursor + 1
        if start + width >= prime_bits:
            width = prime_bits - start - 1
        if width <= 0:
            raise ValueError(f"scaled segment {name} does not fit")
        out[name] = Segment(start, width, name)
        cursor = start + width
    return out


def shape_variables(shape: str, segments: dict[str, Segment], prime_bits: int) -> tuple[tuple[Segment, ...], tuple[Segment, ...]]:
    low_a = segments["low_a"]
    low_b = segments["low_b"]
    high_a = segments["high_a"]
    high_b = segments["high_b"]
    high_c = segments["high_c"]
    edge_low = segments["edge_low"]
    edge_high = segments["edge_high"]
    low_group = Segment(low_a.start, low_b.stop - low_a.start, "low_group")
    high_group = Segment(high_a.start, high_c.stop - high_a.start, "high_group")
    if shape == "S0_grouped_2var":
        return (low_group, high_group), (edge_low, edge_high)
    if shape == "S2_low_exact_high_grouped":
        return (low_a, low_b, high_group), (edge_low, edge_high)
    if shape == "partial_B":
        high_tail = Segment(high_a.start, prime_bits - high_a.start, "high_tail")
        return (low_a, high_tail), (edge_low, low_b)
    if shape == "low_tail_univariate":
        high_tail = Segment(high_a.start, prime_bits - high_a.start, "high_tail")
        return (high_tail,), (edge_low, low_a, low_b)
    raise ValueError(f"unknown shape {shape}")


def make_instance(prime_bits: int, seed: int):
    from sage.all import next_prime

    rng = random.Random(seed)
    primes = []
    while len(primes) < 2:
        candidate = (1 << (prime_bits - 1)) | rng.getrandbits(prime_bits - 1) | 1
        prime = int(next_prime(candidate))
        if prime.bit_length() == prime_bits and prime not in primes:
            primes.append(prime)
    p, q = primes
    return p, q, p * q


def variable_name(seg: Segment) -> str:
    return f"v{seg.start}_{seg.width}"


def build_model(shape: str, prime_bits: int, p: int, n: int):
    segments = scaled_segments(prime_bits)
    variables, fixed_segments = shape_variables(shape, segments, prime_bits)
    base = p
    for seg in variables:
        base &= ~seg.mask
    fixed = []
    for seg in fixed_segments:
        value = (p >> seg.start) & ((1 << seg.width) - 1)
        base &= ~seg.mask
        base |= value << seg.start
        fixed.append((seg, value))
    true_values = []
    for seg in variables:
        value = (p >> seg.start) & ((1 << seg.width) - 1)
        true_values.append((seg, value))
    return {
        "shape": shape,
        "prime_bits": prime_bits,
        "n": n,
        "p": p,
        "base": base,
        "segments": segments,
        "variables": variables,
        "fixed": fixed,
        "true_values": true_values,
        "variable_bits": sum(seg.width for seg in variables),
        "fixed_bits": sum(seg.width for seg, _ in fixed),
    }


def model_record(model, status: str, roots_count: int | None = None, elapsed: float | None = None):
    record = {
        "status": status,
        "shape": model["shape"],
        "prime_bits": model["prime_bits"],
        "n_bits": model["n"].bit_length(),
        "fixed_bits": model["fixed_bits"],
        "variable_bits": model["variable_bits"],
        "fixed": [
            {"name": seg.name, "start": seg.start, "width": seg.width, "value": value}
            for seg, value in model["fixed"]
        ],
        "variables": [
            {"name": variable_name(seg), "segment": seg.name, "start": seg.start, "width": seg.width}
            for seg in model["variables"]
        ],
    }
    if roots_count is not None:
        record["roots"] = roots_count
    if elapsed is not None:
        record["elapsed"] = round(elapsed, 3)
    return record


def print_record(record, as_json: bool) -> None:
    if as_json:
        print(json.dumps(record, sort_keys=True), flush=True)
        return
    variables = ", ".join(f"{item['segment']}:{item['width']}" for item in record["variables"])
    fixed = ", ".join(f"{item['name']}:{item['width']}={item['value']:#x}" for item in record["fixed"])
    parts = [
        f"status={record['status']}",
        f"shape={record['shape']}",
        f"prime_bits={record['prime_bits']}",
        f"variable_bits={record['variable_bits']}",
        f"fixed_bits={record['fixed_bits']}",
        f"fixed=[{fixed}]",
        f"variables=[{variables}]",
    ]
    if "roots" in record:
        parts.append(f"roots={record['roots']}")
    if "elapsed" in record:
        parts.append(f"elapsed={record['elapsed']}")
    print(" ".join(parts), flush=True)


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["dry-run", "cuso", "self-test", "list-shapes"], default="dry-run")
    ap.add_argument("--shape", choices=["all", *SHAPES], default="all")
    ap.add_argument("--prime-bits", type=int, default=128)
    ap.add_argument("--seed", type=int, default=20260705)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--cuso-log", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--cuso-graph", choices=["auto", "on", "off"], default="off")
    ap.add_argument("--cuso-no-intermediate", action="store_true")
    ap.add_argument("--cuso-allow-partial", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    shapes = SHAPES if args.shape == "all" else (args.shape,)

    if args.mode == "self-test":
        p = (1 << (args.prime_bits - 1)) | 0x123456789abcdef
        n = p * ((1 << (args.prime_bits - 1)) | 0xfedcba987654321)
        for shape in SHAPES:
            model = build_model(shape, args.prime_bits, p, n)
            rebuilt = model["base"]
            for seg, value in model["true_values"]:
                rebuilt += value << seg.start
            assert rebuilt == p, (shape, hex(rebuilt), hex(p))
            assert model["variable_bits"] > 0, shape
            assert model["fixed_bits"] > 0, shape
        print("self-test ok")
        return 0

    if args.mode == "list-shapes":
        p = (1 << (args.prime_bits - 1)) | 0x123456789abcdef
        n = p * ((1 << (args.prime_bits - 1)) | 0xfedcba987654321)
        for shape in shapes:
            print_record(model_record(build_model(shape, args.prime_bits, p, n), "dry_run"), args.json)
        return 0

    p, q, n = make_instance(args.prime_bits, args.seed)
    if args.cuso_log:
        logging.basicConfig(
            level=getattr(logging, args.cuso_log),
            format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        )

    exit_status = 0 if args.mode == "dry-run" else 1
    for shape in shapes:
        model = build_model(shape, args.prime_bits, p, n)
        if args.mode == "dry-run":
            print_record(model_record(model, "dry_run"), args.json)
            continue

        from sage.all import PolynomialRing, ZZ
        import cuso

        _patch_cuso_boundset_for_sage109()
        names = [variable_name(seg) for seg in model["variables"]]
        R = PolynomialRing(ZZ, names)
        gens = R.gens()
        f = R(model["base"])
        bounds = {}
        for gen, seg in zip(gens, model["variables"]):
            f += R(1 << seg.start) * gen
            bounds[gen] = (0, 1 << seg.width)

        graph = None if args.cuso_graph == "auto" else args.cuso_graph == "on"
        t0 = time.time()
        roots = cuso_find_small_roots_compat(
            cuso,
            [f],
            bounds,
            modulus="p",
            modulus_multiple=n,
            modulus_lower_bound=1 << (args.prime_bits - 1),
            modulus_upper_bound=1 << args.prime_bits,
            use_graph_optimization=graph,
            use_intermediate_sizes=not args.cuso_no_intermediate,
            allow_partial_solutions=args.cuso_allow_partial,
        )
        elapsed = time.time() - t0
        status = "missed_true_root"
        for root in roots:
            candidate = model["base"]
            if isinstance(root, dict) and "p" in root:
                candidate = int(root["p"])
            else:
                ok = True
                for gen, seg in zip(gens, model["variables"]):
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
            if 1 < candidate < n and n % candidate == 0:
                status = "factor"
                break
            status = "candidate"
        print_record(model_record(model, status, len(roots), elapsed), args.json)
        if status == "factor":
            exit_status = 0
        else:
            exit_status = 1
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
