#!/usr/bin/env python3
"""Mixed-shape cuso runner for ct07_cuso_mixed_shape_search.

The edge 4-bit blocks are fixed with the same candidate-id convention as
`src/solve7_main.py`: `cid = 16*high_edge + low_edge`.
The remaining unknown region is represented by grouped, exact, or hybrid cuso
variables so that the same candidate can be compared across shape choices.
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

from solve7_main import (  # noqa: E402
    LEAK,
    MASK,
    N,
    _patch_cuso_boundset_for_sage109,
    cuso_find_small_roots_compat,
    testp,
)


@dataclass(frozen=True, order=True)
class Segment:
    start: int
    width: int

    @property
    def stop(self) -> int:
        return self.start + self.width

    @property
    def mask(self) -> int:
        return ((1 << self.width) - 1) << self.start

    def label(self) -> str:
        return f"{self.start}:{self.width}"


EDGE_LOW = Segment(150, 4)
EDGE_HIGH = Segment(920, 4)
UNKNOWN_BLOCKS = (
    EDGE_LOW,
    Segment(265, 84),
    Segment(362, 58),
    Segment(600, 69),
    Segment(682, 87),
    Segment(784, 46),
    EDGE_HIGH,
)
NON_EDGE_UNKNOWN_BLOCKS = tuple(seg for seg in UNKNOWN_BLOCKS if seg not in (EDGE_LOW, EDGE_HIGH))

SHAPES = {
    "S0_grouped_2var": (
        Segment(265, 155),
        Segment(600, 230),
    ),
    "S1_exact_5var": (
        Segment(265, 84),
        Segment(362, 58),
        Segment(600, 69),
        Segment(682, 87),
        Segment(784, 46),
    ),
    "S2_low_exact_high_grouped": (
        Segment(265, 84),
        Segment(362, 58),
        Segment(600, 230),
    ),
    "S3_low_grouped_high_exact": (
        Segment(265, 155),
        Segment(600, 69),
        Segment(682, 87),
        Segment(784, 46),
    ),
    "S4_low_exact_high_mixed": (
        Segment(265, 84),
        Segment(362, 58),
        Segment(600, 169),
        Segment(784, 46),
    ),
}


def variable_name(seg: Segment) -> str:
    return f"v{seg.start}_{seg.width}"


def bits_inside(block: Segment, carrier: Segment) -> int:
    left = max(block.start, carrier.start)
    right = min(block.stop, carrier.stop)
    return max(0, right - left)


def unknown_bits_covered_by(segments: tuple[Segment, ...]) -> int:
    return sum(sum(bits_inside(block, seg) for block in NON_EDGE_UNKNOWN_BLOCKS) for seg in segments)


def decode_edge_id(edge_id: int) -> list[tuple[Segment, int]]:
    if not 0 <= edge_id < 256:
        raise ValueError("edge_id must fit in 8 bits")
    low = edge_id & 15
    high = edge_id >> 4
    return [(EDGE_LOW, low), (EDGE_HIGH, high)]


def build_model(
    shape: str,
    edge_id: int,
    upper_variables: list[str] | None = None,
    upper_low_variables: bool = False,
    upper_all_variables: bool = False,
):
    variables = SHAPES[shape]
    edges = decode_edge_id(edge_id)
    upper_segments = set()
    for raw in upper_variables or []:
        try:
            start_s, width_s = raw.split(":", 1)
            upper_segments.add(Segment(int(start_s, 0), int(width_s, 0)))
        except ValueError as ex:
            raise ValueError(f"invalid upper variable spec {raw!r}; expected START:WIDTH") from ex
    unknown_upper = upper_segments - set(variables)
    if unknown_upper:
        labels = ", ".join(seg.label() for seg in sorted(unknown_upper))
        raise ValueError(f"upper variable spec does not match model variables: {labels}")
    variable_origins = {}
    for seg in variables:
        origin = "lower"
        if upper_all_variables or (upper_low_variables and seg.stop <= 600) or seg in upper_segments:
            origin = "upper"
        variable_origins[seg] = origin

    base = LEAK
    for seg in variables:
        base &= ~seg.mask
    for seg, value in edges:
        base &= ~seg.mask
        base |= value << seg.start

    variable_bits = sum(seg.width for seg in variables)
    unknown_bits = unknown_bits_covered_by(variables)
    absorbed_known_bits = variable_bits - unknown_bits
    return {
        "shape": shape,
        "edge_id": edge_id,
        "fixed": edges,
        "variables": variables,
        "variable_origins": variable_origins,
        "base": base,
        "fixed_bits": sum(seg.width for seg, _ in edges),
        "variable_bits": variable_bits,
        "unknown_variable_bits": unknown_bits,
        "absorbed_known_bits": absorbed_known_bits,
    }


def model_record(model, status: str, roots_count: int | None = None, elapsed: float | None = None):
    record = {
        "status": status,
        "shape": model["shape"],
        "edge_id": model["edge_id"],
        "fixed_bits": model["fixed_bits"],
        "variable_bits": model["variable_bits"],
        "unknown_variable_bits": model["unknown_variable_bits"],
        "absorbed_known_bits": model["absorbed_known_bits"],
        "fixed": [
            {"start": seg.start, "width": seg.width, "value": value}
            for seg, value in model["fixed"]
        ],
        "variables": [
            {
                "name": variable_name(seg),
                "start": seg.start,
                "width": seg.width,
                "origin": model["variable_origins"][seg],
            }
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
    fixed = ", ".join(f"{item['start']}:{item['width']}={item['value']:#x}" for item in record["fixed"])
    variables = ", ".join(f"{item['name']}:{item['width']}:{item['origin']}" for item in record["variables"])
    parts = [
        f"status={record['status']}",
        f"shape={record['shape']}",
        f"edge_id={record['edge_id']}",
        f"fixed_bits={record['fixed_bits']}",
        f"variable_bits={record['variable_bits']}",
        f"unknown_variable_bits={record['unknown_variable_bits']}",
        f"absorbed_known_bits={record['absorbed_known_bits']}",
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
    ap.add_argument("--shape", choices=sorted(SHAPES), default="S2_low_exact_high_grouped")
    ap.add_argument("--a", type=int, default=0, help="start edge-id, inclusive")
    ap.add_argument("--b", type=int, help="stop edge-id, exclusive; defaults to a+1")
    ap.add_argument("--upper-variable", action="append", default=[], help="model variable START:WIDTH to use upper origin")
    ap.add_argument("--upper-low-variables", action="store_true", help="use upper origin for variables ending before bit 600")
    ap.add_argument("--upper-all-variables", action="store_true", help="use upper origin for every model variable")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--cuso-log", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--cuso-graph", choices=["auto", "on", "off"], default="auto")
    ap.add_argument("--cuso-no-intermediate", action="store_true")
    ap.add_argument("--cuso-allow-partial", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    if args.mode == "self-test":
        expected = {
            "S0_grouped_2var": (385, 344, 41),
            "S1_exact_5var": (344, 344, 0),
            "S2_low_exact_high_grouped": (372, 344, 28),
            "S3_low_grouped_high_exact": (357, 344, 13),
            "S4_low_exact_high_mixed": (357, 344, 13),
        }
        for shape, values in expected.items():
            model = build_model(shape, 0)
            got = (
                model["variable_bits"],
                model["unknown_variable_bits"],
                model["absorbed_known_bits"],
            )
            assert got == values, (shape, got, values)
            assert model["fixed_bits"] == 8, (shape, model["fixed_bits"])
        upper = build_model("S2_low_exact_high_grouped", 0, upper_low_variables=True)
        assert [upper["variable_origins"][seg] for seg in upper["variables"]] == ["upper", "upper", "lower"]
        explicit = build_model("S0_grouped_2var", 0, upper_variables=["265:155"])
        assert explicit["variable_origins"][Segment(265, 155)] == "upper"
        print("self-test ok")
        return 0

    if args.mode == "list-shapes":
        for shape in sorted(SHAPES):
            model = build_model(shape, 0, args.upper_variable, args.upper_low_variables, args.upper_all_variables)
            variables = ", ".join(seg.label() for seg in model["variables"])
            print(
                f"{shape}: variables={variables} "
                f"variable_bits={model['variable_bits']} "
                f"unknown_variable_bits={model['unknown_variable_bits']} "
                f"absorbed_known_bits={model['absorbed_known_bits']}"
            )
        return 0

    stop = args.b if args.b is not None else args.a + 1
    start = max(0, args.a)
    stop = min(256, stop)
    if start >= stop:
        print(f"empty edge-id range {start}:{stop}", file=sys.stderr)
        return 2

    if args.cuso_log:
        logging.basicConfig(
            level=getattr(logging, args.cuso_log),
            format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        )

    exit_status = 1
    for edge_id in range(start, stop):
        model = build_model(
            args.shape,
            edge_id,
            args.upper_variable,
            args.upper_low_variables,
            args.upper_all_variables,
        )
        if args.mode == "dry-run":
            print_record(model_record(model, "dry_run"), args.json)
            exit_status = 0
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
            coeff = R(1 << seg.start)
            if model["variable_origins"][seg] == "upper":
                f += coeff * ((1 << seg.width) - 1)
                f -= coeff * gen
            else:
                f += coeff * gen
            bounds[gen] = (0, 1 << seg.width)

        graph = None if args.cuso_graph == "auto" else args.cuso_graph == "on"
        t0 = time.time()
        roots = cuso_find_small_roots_compat(
            cuso,
            [f],
            bounds,
            modulus="p",
            modulus_multiple=N,
            modulus_lower_bound=1 << 1023,
            modulus_upper_bound=1 << 1024,
            use_graph_optimization=graph,
            use_intermediate_sizes=not args.cuso_no_intermediate,
            allow_partial_solutions=args.cuso_allow_partial,
        )
        elapsed = time.time() - t0
        status = "soft_no_root"
        if roots:
            status = "candidate"
        for root in roots:
            if isinstance(root, dict) and "p" in root:
                p = int(root["p"])
            else:
                p = model["base"]
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
                    if model["variable_origins"][seg] == "upper":
                        value = ((1 << seg.width) - 1) - value
                    p += value << seg.start
                if not ok:
                    continue
            if testp(p):
                status = "factor"
                break
        print_record(model_record(model, status, len(roots), elapsed), args.json)
        if status == "factor":
            return 0
        exit_status = 1
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
