#!/usr/bin/env python3
"""Partial-low600 cuso oracle for broad-clause experiments.

This is the first implementation target for ct07_partial_low600_cuso_broad_clause.
It builds unknown-divisor cuso relations of the form

    p = C + sum(2^start * u_start) + 2^600 * z600

where only part of the low600 unknown bits are fixed.  A no-root result is
reported as soft_no_root by default; do not turn it into a SAT hard clause until
an external soundness gate has validated the shape.
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


LOW600_UNKNOWN_BLOCKS = (
    Segment(150, 4),
    Segment(265, 84),
    Segment(362, 58),
)
HIGH_TAIL = Segment(600, 424)

SHAPES = {
    "none": (),
    "A": (Segment(362, 58),),
    "B": (Segment(150, 4), Segment(362, 58)),
    "C": (Segment(265, 84),),
    "D_265_64_362_8": (Segment(265, 64), Segment(362, 8)),
    "D_265_48_362_16": (Segment(265, 48), Segment(362, 16)),
    "E_150_4_265_84_362_16": (Segment(150, 4), Segment(265, 84), Segment(362, 16)),
    "E_150_4_265_84_362_32": (Segment(150, 4), Segment(265, 84), Segment(362, 32)),
    "F_150_4_265_84_362_50": (Segment(150, 4), Segment(265, 84), Segment(362, 50)),
    "F_150_4_265_84_362_54": (Segment(150, 4), Segment(265, 84), Segment(362, 54)),
    "F_150_4_265_84_362_58": (Segment(150, 4), Segment(265, 84), Segment(362, 58)),
}


def parse_fixed(text: str) -> tuple[Segment, int]:
    try:
        left, value_s = text.split("=", 1)
        start_s, width_s = left.split(":", 1)
        seg = Segment(int(start_s, 0), int(width_s, 0))
        value = int(value_s, 0)
    except ValueError as ex:
        raise argparse.ArgumentTypeError(
            "expected START:WIDTH=VALUE, for example 150:4=0xa"
        ) from ex
    if seg.width <= 0:
        raise argparse.ArgumentTypeError("fixed segment width must be positive")
    if not 0 <= value < (1 << seg.width):
        raise argparse.ArgumentTypeError(f"value does not fit in {seg.width} bits")
    return seg, value


def segment_inside_any(seg: Segment, blocks: tuple[Segment, ...]) -> bool:
    return any(block.start <= seg.start and seg.stop <= block.stop for block in blocks)


def validate_fixed_segments(fixed: list[tuple[Segment, int]]) -> None:
    ordered = sorted(fixed, key=lambda item: item[0].start)
    for seg, value in ordered:
        if not segment_inside_any(seg, LOW600_UNKNOWN_BLOCKS):
            raise ValueError(f"fixed segment {seg.label()} is not inside low600 unknown blocks")
        if not 0 <= value < (1 << seg.width):
            raise ValueError(f"value for {seg.label()} does not fit")
    for (left, _), (right, _) in zip(ordered, ordered[1:]):
        if right.start < left.stop:
            raise ValueError(f"fixed segments overlap: {left.label()} and {right.label()}")


def decode_fixed_id(shape_segments: tuple[Segment, ...], fixed_id: int) -> list[tuple[Segment, int]]:
    if fixed_id < 0:
        raise ValueError("fixed_id must be nonnegative")
    fixed = []
    tmp = fixed_id
    for seg in shape_segments:
        value = tmp & ((1 << seg.width) - 1)
        tmp >>= seg.width
        fixed.append((seg, value))
    if tmp:
        raise ValueError("fixed_id exceeds selected shape space")
    return fixed


def variable_low_segments(fixed: list[tuple[Segment, int]]) -> list[Segment]:
    fixed_segments = sorted((seg for seg, _ in fixed), key=lambda seg: seg.start)
    variables: list[Segment] = []
    for block in LOW600_UNKNOWN_BLOCKS:
        cursor = block.start
        inside = [
            seg
            for seg in fixed_segments
            if block.start <= seg.start and seg.stop <= block.stop
        ]
        for seg in inside:
            if cursor < seg.start:
                variables.append(Segment(cursor, seg.start - cursor))
            cursor = seg.stop
        if cursor < block.stop:
            variables.append(Segment(cursor, block.stop - cursor))
    return variables


def build_model(
    shape: str,
    fixed_id: int,
    extra_fixed: list[tuple[Segment, int]],
    upper_variables: list[str] | None = None,
    upper_low_variables: bool = False,
    upper_all_variables: bool = False,
):
    fixed = decode_fixed_id(SHAPES[shape], fixed_id) + list(extra_fixed)
    validate_fixed_segments(fixed)
    low_variables = variable_low_segments(fixed)
    variables = low_variables + [HIGH_TAIL]
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
        if upper_all_variables or (upper_low_variables and seg != HIGH_TAIL) or seg in upper_segments:
            origin = "upper"
        variable_origins[seg] = origin

    base = LEAK
    for block in LOW600_UNKNOWN_BLOCKS:
        base &= ~block.mask
    base &= ~HIGH_TAIL.mask
    for seg, value in fixed:
        base |= value << seg.start

    fixed_bits = sum(seg.width for seg, _ in fixed)
    variable_bits = sum(seg.width for seg in variables)
    fixed_space_bits = sum(seg.width for seg in SHAPES[shape])
    return {
        "shape": shape,
        "fixed_id": fixed_id,
        "fixed": fixed,
        "variables": variables,
        "variable_origins": variable_origins,
        "base": base,
        "fixed_bits": fixed_bits,
        "variable_bits": variable_bits,
        "fixed_space_bits": fixed_space_bits,
    }


def shape_space(shape: str) -> int:
    return 1 << sum(seg.width for seg in SHAPES[shape])


def model_record(model, status: str, roots_count: int | None = None, elapsed: float | None = None):
    record = {
        "status": status,
        "shape": model["shape"],
        "fixed_id": model["fixed_id"],
        "fixed_bits": model["fixed_bits"],
        "variable_bits": model["variable_bits"],
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


def variable_name(seg: Segment) -> str:
    if seg == HIGH_TAIL:
        return "z600"
    return f"u{seg.start}_{seg.width}"


def print_record(record, as_json: bool) -> None:
    if as_json:
        print(json.dumps(record, sort_keys=True), flush=True)
        return
    fixed = ", ".join(f"{item['start']}:{item['width']}={item['value']:#x}" for item in record["fixed"])
    variables = ", ".join(f"{item['name']}:{item['width']}:{item['origin']}" for item in record["variables"])
    parts = [
        f"status={record['status']}",
        f"shape={record['shape']}",
        f"fixed_id={record['fixed_id']}",
        f"fixed_bits={record['fixed_bits']}",
        f"variable_bits={record['variable_bits']}",
        f"fixed=[{fixed}]",
        f"variables=[{variables}]",
    ]
    if "roots" in record:
        parts.append(f"roots={record['roots']}")
    if "elapsed" in record:
        parts.append(f"elapsed={record['elapsed']}")
    print(" ".join(parts), flush=True)


def build_sage_relation(model):
    from sage.all import PolynomialRing, ZZ

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
    return f, bounds, gens


def root_value(root, gen):
    for key in (gen, str(gen)):
        try:
            if key in root:
                return int(root[key])
        except Exception:
            pass
    return None


def reconstruct_p(model, root, gens) -> int | None:
    if isinstance(root, dict) and "p" in root:
        return int(root["p"])
    p = model["base"]
    for gen, seg in zip(gens, model["variables"]):
        value = root_value(root, gen)
        if value is None:
            return None
        if model["variable_origins"][seg] == "upper":
            value = ((1 << seg.width) - 1) - value
        p += value << seg.start
    return p


def run_cuso(model, args) -> tuple[str, int, float]:
    if args.cuso_log:
        logging.basicConfig(
            level=getattr(logging, args.cuso_log),
            format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        )
    import cuso

    _patch_cuso_boundset_for_sage109()
    f, bounds, gens = build_sage_relation(model)
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
    for root in roots:
        p = reconstruct_p(model, root, gens)
        if p is not None and testp(p):
            return "factor", len(roots), elapsed
    if roots:
        return "candidate", len(roots), elapsed
    return "soft_no_root", 0, elapsed


def run_leaf_certify(model, args) -> tuple[list[dict], dict]:
    low_variables = [seg for seg in model["variables"] if seg != HIGH_TAIL]
    low_bits = sum(seg.width for seg in low_variables)
    completions = 1 << low_bits
    leaf_start = max(0, getattr(args, "leaf_start", 0) or 0)
    leaf_stop_arg = getattr(args, "leaf_stop", None)
    leaf_stop = completions if leaf_stop_arg is None else min(completions, leaf_stop_arg)
    if leaf_stop < leaf_start:
        leaf_stop = leaf_start
    selected_completions = leaf_stop - leaf_start
    max_completions = getattr(args, "leaf_max_completions", 16)
    if selected_completions > max_completions:
        return [], {
            **model_record(model, "leaf_too_wide"),
            "leaf_low_variable_bits": low_bits,
            "leaf_completions": completions,
            "leaf_range_start": leaf_start,
            "leaf_range_stop": leaf_stop,
            "leaf_range_completions": selected_completions,
            "leaf_max_completions": max_completions,
        }
    if selected_completions == 0:
        return [], {
            **model_record(model, "leaf_empty_range"),
            "leaf_low_variable_bits": low_bits,
            "leaf_completions": completions,
            "leaf_range_start": leaf_start,
            "leaf_range_stop": leaf_stop,
            "leaf_range_completions": selected_completions,
            "leaf_max_completions": max_completions,
        }

    status_counts: dict[str, int] = {}
    elapsed_total = 0.0
    leaf_records = []
    high_tail_upper = model["variable_origins"].get(HIGH_TAIL) == "upper"
    for completion_id in range(leaf_start, leaf_stop):
        tmp = completion_id
        extra = []
        for seg in low_variables:
            value = tmp & ((1 << seg.width) - 1)
            tmp >>= seg.width
            extra.append((seg, value))
        leaf_model = build_model(
            "none",
            0,
            model["fixed"] + extra,
            [],
            False,
            high_tail_upper,
        )
        status, roots_count, elapsed = run_cuso(leaf_model, args)
        elapsed_total += elapsed
        status_counts[status] = status_counts.get(status, 0) + 1
        leaf_record = {
            **model_record(leaf_model, status, roots_count, elapsed),
            "parent_shape": model["shape"],
            "parent_fixed_id": model["fixed_id"],
            "leaf_completion_id": completion_id,
            "leaf_completion_count": completions,
            "leaf_range_start": leaf_start,
            "leaf_range_stop": leaf_stop,
        }
        leaf_records.append(leaf_record)
        if status == "factor":
            return leaf_records, {
                **model_record(model, "factor", elapsed=elapsed_total),
                "leaf_low_variable_bits": low_bits,
                "leaf_completions_checked": len(leaf_records),
                "leaf_completions": completions,
                "leaf_range_start": leaf_start,
                "leaf_range_stop": leaf_stop,
                "leaf_range_completions": selected_completions,
                "leaf_status_counts": status_counts,
            }
        if status != "soft_no_root":
            return leaf_records, {
                **model_record(model, "leaf_inconclusive", elapsed=elapsed_total),
                "leaf_low_variable_bits": low_bits,
                "leaf_completions_checked": len(leaf_records),
                "leaf_completions": completions,
                "leaf_range_start": leaf_start,
                "leaf_range_stop": leaf_stop,
                "leaf_range_completions": selected_completions,
                "leaf_status_counts": status_counts,
            }

    status = "leaf_exhausted_soft_no_root"
    if leaf_start != 0 or leaf_stop != completions:
        status = "leaf_range_exhausted_soft_no_root"
    return leaf_records, {
        **model_record(model, status, elapsed=elapsed_total),
        "leaf_low_variable_bits": low_bits,
        "leaf_completions_checked": len(leaf_records),
        "leaf_completions": completions,
        "leaf_range_start": leaf_start,
        "leaf_range_stop": leaf_stop,
        "leaf_range_completions": selected_completions,
        "leaf_status_counts": status_counts,
    }


def run_self_test() -> None:
    expected = {
        "A": (58, 512),
        "B": (62, 508),
        "C": (84, 486),
        "D_265_64_362_8": (72, 498),
        "D_265_48_362_16": (64, 506),
        "E_150_4_265_84_362_16": (104, 466),
        "E_150_4_265_84_362_32": (120, 450),
        "F_150_4_265_84_362_50": (138, 432),
        "F_150_4_265_84_362_54": (142, 428),
        "F_150_4_265_84_362_58": (146, 424),
    }
    for shape, (fixed_bits, variable_bits) in expected.items():
        model = build_model(shape, 0, [])
        assert model["fixed_bits"] == fixed_bits, (shape, model["fixed_bits"])
        assert model["variable_bits"] == variable_bits, (shape, model["variable_bits"])
    custom = build_model("none", 0, [parse_fixed("150:4=0xf"), parse_fixed("265:8=0xaa")])
    assert custom["fixed_bits"] == 12
    assert custom["variable_bits"] == 76 + 58 + 424
    upper = build_model("B", 0, [], upper_low_variables=True)
    assert [upper["variable_origins"][seg] for seg in upper["variables"]] == ["upper", "lower"]
    explicit = build_model("B", 0, [], upper_variables=["265:84"])
    assert explicit["variable_origins"][Segment(265, 84)] == "upper"
    f54 = build_model("F_150_4_265_84_362_54", 0, [])
    assert [seg.label() for seg in f54["variables"]] == ["416:4", "600:424"]
    f50 = build_model("F_150_4_265_84_362_50", 0, [])
    assert [seg.label() for seg in f50["variables"]] == ["412:8", "600:424"]
    range_args = argparse.Namespace(
        cuso_log=None,
        cuso_graph="off",
        cuso_no_intermediate=False,
        cuso_allow_partial=False,
        leaf_max_completions=16,
        leaf_start=0,
        leaf_stop=17,
    )
    _, too_wide = run_leaf_certify(f50, range_args)
    assert too_wide["status"] == "leaf_too_wide"
    range_args.leaf_start = 2
    range_args.leaf_stop = 16
    _, empty_or_range = run_leaf_certify(build_model("F_150_4_265_84_362_58", 0, []), range_args)
    assert empty_or_range["status"] == "leaf_empty_range"
    print("self-test ok")


def list_shapes() -> None:
    for shape in sorted(SHAPES):
        model = build_model(shape, 0, [])
        fixed = ", ".join(seg.label() for seg in SHAPES[shape]) or "(none)"
        variables = ", ".join(seg.label() for seg in model["variables"])
        print(
            f"{shape}: fixed={fixed} fixed_bits={model['fixed_bits']} "
            f"variable_bits={model['variable_bits']} variables={variables}"
        )


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--mode",
        choices=["dry-run", "cuso", "leaf-certify", "self-test", "list-shapes"],
        default="dry-run",
    )
    ap.add_argument("--shape", choices=sorted(SHAPES), default="B")
    ap.add_argument("--a", type=int, default=0, help="start fixed-id, inclusive")
    ap.add_argument("--b", type=int, help="stop fixed-id, exclusive; defaults to a+1")
    ap.add_argument("--fixed", action="append", type=parse_fixed, default=[], help="extra fixed START:WIDTH=VALUE")
    ap.add_argument("--upper-variable", action="append", default=[], help="model variable START:WIDTH to use upper origin")
    ap.add_argument("--upper-low-variables", action="store_true", help="use upper origin for low600 variables only")
    ap.add_argument("--upper-all-variables", action="store_true", help="use upper origin for every variable, including z600")
    ap.add_argument("--results-jsonl", help="append leaf-certify records to this JSONL file")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--cuso-log", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--cuso-graph", choices=["auto", "on", "off"], default="auto")
    ap.add_argument("--cuso-no-intermediate", action="store_true")
    ap.add_argument("--cuso-allow-partial", action="store_true")
    ap.add_argument(
        "--leaf-max-completions",
        type=int,
        default=16,
        help="maximum low-variable completions to enumerate in leaf-certify mode",
    )
    ap.add_argument("--leaf-start", type=int, default=0, help="leaf completion start, inclusive")
    ap.add_argument("--leaf-stop", type=int, help="leaf completion stop, exclusive")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.mode == "self-test":
        run_self_test()
        return 0
    if args.mode == "list-shapes":
        list_shapes()
        return 0

    stop = args.b if args.b is not None else args.a + 1
    space = shape_space(args.shape)
    start = max(0, args.a)
    stop = min(space, stop)
    if start >= stop:
        print(f"empty fixed-id range {start}:{stop} for shape {args.shape}", file=sys.stderr)
        return 2

    exit_status = 1
    for fixed_id in range(start, stop):
        model = build_model(
            args.shape,
            fixed_id,
            args.fixed,
            args.upper_variable,
            args.upper_low_variables,
            args.upper_all_variables,
        )
        if args.mode == "dry-run":
            print_record(model_record(model, "dry_run"), args.json)
            exit_status = 0
            continue
        if args.mode == "leaf-certify":
            leaf_records, summary = run_leaf_certify(model, args)
            for leaf_record in leaf_records:
                print_record(leaf_record, args.json)
                if args.results_jsonl:
                    path = Path(args.results_jsonl)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with path.open("a", encoding="utf-8") as out:
                        out.write(json.dumps(leaf_record, sort_keys=True) + "\n")
            print_record(summary, args.json)
            if args.results_jsonl:
                path = Path(args.results_jsonl)
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("a", encoding="utf-8") as out:
                    out.write(json.dumps(summary, sort_keys=True) + "\n")
            if summary["status"] == "factor":
                return 0
            if summary["status"] == "leaf_too_wide":
                return 2
            exit_status = 1
            continue
        status, roots_count, elapsed = run_cuso(model, args)
        print_record(model_record(model, status, roots_count, elapsed), args.json)
        if status == "factor":
            return 0
        if status in ("candidate", "soft_no_root"):
            exit_status = 1
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
