#!/usr/bin/env python3
"""Focus-group style HM lattice profiler for planted low-free models.

This is not a solver.  It builds a small Herrmann-May style lattice for the
downscaled planted low-free model and asks which reduced rows vanish at the
known true root, while recording the input shift rows used by LLL's
transformation matrix.  The goal is to identify useful shift families before
trying row-pruned fallback lattices on the full challenge shape.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ASSET_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("DOT_SAGE", str(ASSET_ROOT / ".sage-tmp"))
SRC_ROOT = ASSET_ROOT / "src"
EXPERIMENT_ROOT = ASSET_ROOT / "experiments"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
if str(EXPERIMENT_ROOT) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_ROOT))

from fpylll import IntegerMatrix, LLL  # noqa: E402
from planted_lowfree_sweep import build_model, make_instance, synthetic_instance  # noqa: E402


class NoVanishingRow(Exception):
    def __init__(self, summary: dict):
        super().__init__("source lattice has no vanishing row")
        self.summary = summary


def centered(value: int, modulus: int) -> int:
    value %= modulus
    return value - modulus if value > modulus // 2 else value


def poly_mul(left: dict[tuple[int, int], int], right: dict[tuple[int, int], int]) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {}
    for (a0, a1), left_coeff in left.items():
        for (b0, b1), right_coeff in right.items():
            mon = (a0 + b0, a1 + b1)
            out[mon] = out.get(mon, 0) + left_coeff * right_coeff
    return {mon: coeff for mon, coeff in out.items() if coeff}


def poly_shift_xy(
    poly: dict[tuple[int, int], int],
    x_amount: int,
    y_amount: int,
    scale: int,
) -> dict[tuple[int, int], int]:
    return {
        (x_degree + x_amount, y_degree + y_amount): coeff * scale
        for (x_degree, y_degree), coeff in poly.items()
        if coeff
    }


def poly_reduce_centered(poly: dict[tuple[int, int], int], modulus: int) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {}
    for mon, coeff in poly.items():
        reduced = centered(coeff, modulus)
        if reduced:
            out[mon] = reduced
    return out


def eval_poly(poly: dict[tuple[int, int], int], x_value: int, y_value: int) -> int:
    x_pows = {0: 1}
    y_pows = {0: 1}
    for x_degree, y_degree in poly:
        if x_degree not in x_pows:
            x_pows[x_degree] = x_value**x_degree
        if y_degree not in y_pows:
            y_pows[y_degree] = y_value**y_degree
    total = 0
    for (x_degree, y_degree), coeff in poly.items():
        total += coeff * x_pows[x_degree] * y_pows[y_degree]
    return total


def family_key(item: dict) -> str:
    if item.get("x_shift", 0):
        return f"{item['power']}:{item['x_shift']}:{item['y_shift']}:{item['n_scale_power']}"
    return f"{item['power']}:{item['y_shift']}:{item['n_scale_power']}"


def build_basis(
    model,
    m: int,
    t: int,
    keep_families: set[str] | None = None,
    construction: str = "y-only",
    x_shift_limit: int = 1,
    x_origin: str = "zero",
):
    if len(model.variables) != 2:
        raise ValueError("focus-group HM profiler expects exactly two variables")
    x_seg, y_seg = model.variables
    x_bound = 1 << x_seg.width
    y_bound = 1 << y_seg.width
    x_coeff = 1 << x_seg.start
    y_coeff = 1 << y_seg.start
    inv_x_coeff = pow(x_coeff, -1, model.n)
    if x_origin == "zero":
        y_term = centered(y_coeff * inv_x_coeff, model.n)
        constant_term = centered(model.base * inv_x_coeff, model.n)
    elif x_origin == "upper":
        y_term = centered(-y_coeff * inv_x_coeff, model.n)
        constant_term = centered(-(model.base + x_coeff * (x_bound - 1)) * inv_x_coeff, model.n)
    else:
        raise ValueError(f"unknown x origin: {x_origin}")
    f = {
        (1, 0): 1,
        (0, 1): y_term,
        (0, 0): constant_term,
    }

    n_pows = [1]
    for _ in range(m):
        n_pows.append(n_pows[-1] * model.n)

    powers = [{(0, 0): 1}]
    current = {(0, 0): 1}
    for power in range(1, m + 1):
        current = poly_reduce_centered(poly_mul(current, f), n_pows[power])
        powers.append(current)

    shifts = []
    for power in range(m + 1):
        n_scale = n_pows[max(t - power, 0)]
        max_shift = m - power
        if construction == "y-only":
            shift_pairs = [(0, y_shift) for y_shift in range(max_shift + 1)]
        elif construction == "x-limited":
            capped_x_shift = min(max(x_shift_limit, 0), max_shift)
            shift_pairs = [
                (x_shift, y_shift)
                for x_shift in range(capped_x_shift + 1)
                for y_shift in range(max_shift - x_shift + 1)
            ]
        elif construction == "total-degree":
            shift_pairs = [
                (x_shift, y_shift)
                for x_shift in range(max_shift + 1)
                for y_shift in range(max_shift - x_shift + 1)
            ]
        else:
            raise ValueError(f"unknown construction: {construction}")
        for x_shift, y_shift in shift_pairs:
            shifts.append(
                {
                    "power": power,
                    "x_shift": x_shift,
                    "y_shift": y_shift,
                    "n_scale_power": max(t - power, 0),
                    "poly": poly_shift_xy(powers[power], x_shift, y_shift, n_scale),
                }
            )
    if keep_families is not None:
        shifts = [item for item in shifts if family_key(item) in keep_families]
        if not shifts:
            raise ValueError("row-pruned basis has no retained shifts")

    monomials = sorted(
        set().union(*(set(item["poly"].keys()) for item in shifts)),
        key=lambda mon: (mon[0] + mon[1], mon[0], mon[1]),
    )
    scales = [(x_bound**x_degree) * (y_bound**y_degree) for x_degree, y_degree in monomials]
    matrix = IntegerMatrix(len(shifts), len(monomials))
    for row_index, item in enumerate(shifts):
        for col_index, mon in enumerate(monomials):
            matrix[row_index, col_index] = item["poly"].get(mon, 0) * scales[col_index]
    return matrix, monomials, scales, shifts


def row_to_poly(matrix: IntegerMatrix, row_index: int, monomials, scales):
    poly: dict[tuple[int, int], int] = {}
    for col_index, mon in enumerate(monomials):
        value = int(matrix[row_index, col_index])
        if not value:
            continue
        scale = scales[col_index]
        if value % scale:
            return None
        poly[mon] = value // scale
    return poly


def row_norm_bits(matrix: IntegerMatrix, row_index: int) -> float:
    norm2 = 0
    for col_index in range(matrix.ncols):
        value = int(matrix[row_index, col_index])
        norm2 += value * value
    return 0.0 if norm2 == 0 else norm2.bit_length() / 2


def parse_int_set(text: str) -> list[int]:
    out = []
    for raw_part in text.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_s, stop_s = part.split("-", 1)
            out.extend(range(int(start_s), int(stop_s) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def profile_model(
    model,
    m: int,
    t: int,
    true_values: list[int],
    row_limit: int,
    keep_families: set[str] | None = None,
    construction: str = "y-only",
    x_shift_limit: int = 1,
    x_origin: str = "zero",
) -> dict:
    matrix, monomials, scales, shifts = build_basis(
        model,
        m,
        t,
        keep_families,
        construction,
        x_shift_limit,
        x_origin,
    )
    transform = IntegerMatrix.identity(matrix.nrows)
    start_time = time.time()
    LLL.reduction(matrix, transform, delta=0.99, method="proved", float_type="mpfr", precision=192)
    elapsed = time.time() - start_time

    rows = []
    limit = matrix.nrows if row_limit <= 0 else min(row_limit, matrix.nrows)
    for row_index in range(limit):
        poly = row_to_poly(matrix, row_index, monomials, scales)
        if poly is None:
            rows.append({"row": row_index, "status": "not_divisible"})
            continue
        if not poly:
            rows.append(
                {
                    "row": row_index,
                    "status": "zero_row",
                    "norm_bits": 0.0,
                    "terms": 0,
                    "degree": 0,
                    "vanishes_at_true_root": False,
                    "value_bits_at_true_root": 0,
                    "support_count": 0,
                    "support": [],
                }
            )
            continue
        support = []
        for input_index in range(transform.ncols):
            coeff = int(transform[row_index, input_index])
            if coeff:
                item = shifts[input_index]
                support.append(
                    {
                        "input": input_index,
                        "coeff": coeff,
                        "power": item["power"],
                        "x_shift": item["x_shift"],
                        "y_shift": item["y_shift"],
                        "n_scale_power": item["n_scale_power"],
                        "family": family_key(item),
                    }
                )
        value_at_root = eval_poly(poly, true_values[0], true_values[1])
        rows.append(
            {
                "row": row_index,
                "norm_bits": round(row_norm_bits(matrix, row_index), 3),
                "terms": len(poly),
                "degree": max((sum(mon) for mon in poly), default=0),
                "vanishes_at_true_root": value_at_root == 0,
                "value_bits_at_true_root": abs(value_at_root).bit_length() if value_at_root else 0,
                "support_count": len(support),
                "support": support,
            }
        )

    return {
        "status": "profile",
        "construction": construction,
        "x_shift_limit": x_shift_limit,
        "x_origin": x_origin,
        "m": m,
        "t": t,
        "lll_elapsed": round(elapsed, 3),
        "basis_rows": matrix.nrows,
        "basis_cols": matrix.ncols,
        "retained_families": None if keep_families is None else sorted(keep_families),
        "rows": rows,
    }


def first_vanishing_source(
    model,
    m: int,
    t: int,
    true_values: list[int],
    construction: str = "y-only",
    x_shift_limit: int = 1,
    x_origin: str = "zero",
) -> tuple[dict, dict]:
    full_summary = profile_model(
        model,
        m,
        t,
        true_values,
        0,
        construction=construction,
        x_shift_limit=x_shift_limit,
        x_origin=x_origin,
    )
    full_vanishing = [row for row in full_summary["rows"] if row.get("vanishes_at_true_root")]
    if not full_vanishing:
        raise NoVanishingRow(full_summary)
    return full_summary, full_vanishing[0]


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["self-test", "profile", "sweep", "seed-sweep", "prune", "drop-sweep"],
        default="profile",
    )
    parser.add_argument("--prime-bits", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260705)
    parser.add_argument("--low-free-width", type=int, default=1)
    parser.add_argument(
        "--low-free-offset",
        type=int,
        help="offset from the scaled low_group start; default keeps the suffix",
    )
    parser.add_argument("--m", type=int, default=3)
    parser.add_argument("--t", type=int, default=1)
    parser.add_argument("--m-values", default="3-8")
    parser.add_argument("--t-values", default="1-4")
    parser.add_argument("--seeds", default="20260705")
    parser.add_argument("--rows", type=int, default=12)
    parser.add_argument("--construction", choices=["y-only", "x-limited", "total-degree"], default="y-only")
    parser.add_argument("--x-shift-limit", type=int, default=1)
    parser.add_argument("--x-origin", choices=["zero", "upper"], default="zero")
    parser.add_argument("--template-families", default="")
    parser.add_argument("--min-support-coeff-bits", type=int, default=0)
    parser.add_argument("--drop-axes", default="power,nscale")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.mode == "seed-sweep":
        template_families = {
            family.strip()
            for family in args.template_families.split(",")
            if family.strip()
        }
        seed_rows = []
        family_counts: dict[str, int] = {}
        success_cases = 0
        total_cases = 0
        last_variables = []
        for seed in parse_int_set(args.seeds):
            p, _q, n = make_instance(args.prime_bits, seed)
            model = build_model(
                args.prime_bits,
                p,
                n,
                args.low_free_width,
                args.low_free_offset,
            )
            if len(model.variables) != 2:
                raise SystemExit("choose a positive --low-free-width so the model has two variables")
            true_values = []
            for segment in model.variables:
                true_values.append((p >> segment.start) & ((1 << segment.width) - 1))
            rebuilt = model.base
            for segment, value in zip(model.variables, true_values):
                rebuilt += value << segment.start
            if rebuilt != p:
                raise AssertionError((hex(rebuilt), hex(p)))
            root_values = list(true_values)
            if args.x_origin == "upper":
                root_values[0] = ((1 << model.variables[0].width) - 1) - true_values[0]
            last_variables = [
                {
                    "name": segment.name,
                    "start": segment.start,
                    "width": segment.width,
                    "true_value": value,
                    "root_value": root_value,
                }
                for segment, value, root_value in zip(model.variables, true_values, root_values)
            ]
            seed_variables = last_variables
            for m_value in parse_int_set(args.m_values):
                for t_value in parse_int_set(args.t_values):
                    if t_value > m_value:
                        continue
                    total_cases += 1
                    summary = profile_model(
                        model,
                        m_value,
                        t_value,
                        root_values,
                        0,
                        template_families or None,
                        construction=args.construction,
                        x_shift_limit=args.x_shift_limit,
                        x_origin=args.x_origin,
                    )
                    vanishing = [row for row in summary["rows"] if row.get("vanishes_at_true_root")]
                    first = vanishing[0] if vanishing else None
                    support_families = []
                    if first is not None:
                        success_cases += 1
                        support_families = sorted({item["family"] for item in first["support"]})
                        for family in support_families:
                            family_counts[family] = family_counts.get(family, 0) + 1
                    row = {
                        "seed": seed,
                        "variables": seed_variables,
                        "m": m_value,
                        "t": t_value,
                        "basis_rows": summary["basis_rows"],
                        "basis_cols": summary["basis_cols"],
                        "lll_elapsed": summary["lll_elapsed"],
                        "vanishing_rows": len(vanishing),
                        "first_vanishing_row": None if first is None else first["row"],
                        "first_vanishing_norm_bits": None if first is None else first["norm_bits"],
                        "first_vanishing_support_count": None if first is None else first["support_count"],
                        "first_vanishing_support_families": support_families[:24],
                    }
                    seed_rows.append(row)
                    if not args.json:
                        first_text = "-"
                        if first is not None:
                            first_text = (
                                f"row{first['row']} norm={first['norm_bits']} "
                                f"support={first['support_count']}"
                            )
                        print(
                            f"seed={seed} construction={args.construction} "
                            f"x_origin={args.x_origin} x_shift_limit={args.x_shift_limit} "
                            f"low_free_offset={model.low_free_offset} "
                            f"m={m_value} t={t_value} "
                            f"template={len(template_families)} "
                            f"dim={summary['basis_rows']}x{summary['basis_cols']} "
                            f"lll={summary['lll_elapsed']} vanish={len(vanishing)} first={first_text}",
                            flush=True,
                        )
        frequent_families = sorted(family_counts.items(), key=lambda item: (-item[1], item[0]))
        result = {
            "status": "seed-sweep",
            "construction": args.construction,
            "x_shift_limit": args.x_shift_limit,
            "x_origin": args.x_origin,
            "template_families": sorted(template_families),
            "prime_bits": args.prime_bits,
            "low_free_width": args.low_free_width,
            "low_free_offset": None if not seed_rows else model.low_free_offset,
            "seeds": parse_int_set(args.seeds),
            "m_values": parse_int_set(args.m_values),
            "t_values": parse_int_set(args.t_values),
            "variables_last_seed": last_variables,
            "total_cases": total_cases,
            "success_cases": success_cases,
            "success_rate": None if not total_cases else round(success_cases / total_cases, 4),
            "frequent_support_families": [
                {"family": family, "count": count} for family, count in frequent_families[:32]
            ],
            "results": seed_rows,
        }
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True))
        else:
            print(
                f"summary cases={total_cases} success={success_cases} "
                f"rate={result['success_rate']} top_families={frequent_families[:16]}"
            )
        return 0

    if args.mode == "self-test":
        p, _q, n = synthetic_instance(args.prime_bits)
    else:
        p, _q, n = make_instance(args.prime_bits, args.seed)
    model = build_model(
        args.prime_bits,
        p,
        n,
        args.low_free_width,
        args.low_free_offset,
    )
    if len(model.variables) != 2:
        raise SystemExit("choose a positive --low-free-width so the model has two variables")

    true_values = []
    for segment in model.variables:
        true_values.append((p >> segment.start) & ((1 << segment.width) - 1))
    rebuilt = model.base
    for segment, value in zip(model.variables, true_values):
        rebuilt += value << segment.start
    if rebuilt != p:
        raise AssertionError((hex(rebuilt), hex(p)))
    root_values = list(true_values)
    if args.x_origin == "upper":
        root_values[0] = ((1 << model.variables[0].width) - 1) - true_values[0]
    if args.mode == "self-test":
        matrix, monomials, scales, shifts = build_basis(
            model,
            args.m,
            args.t,
            construction=args.construction,
            x_shift_limit=args.x_shift_limit,
            x_origin=args.x_origin,
        )
        assert matrix.nrows == len(shifts)
        assert matrix.ncols == len(monomials) == len(scales)
        print("self-test ok")
        return 0

    variables = [
        {"name": segment.name, "start": segment.start, "width": segment.width, "true_value": value, "root_value": root_value}
        for segment, value, root_value in zip(model.variables, true_values, root_values)
    ]

    if args.mode == "sweep":
        sweep_rows = []
        for m_value in parse_int_set(args.m_values):
            for t_value in parse_int_set(args.t_values):
                if t_value > m_value:
                    continue
                summary = profile_model(
                    model,
                    m_value,
                    t_value,
                    root_values,
                    0,
                    construction=args.construction,
                    x_shift_limit=args.x_shift_limit,
                    x_origin=args.x_origin,
                )
                vanishing = [row for row in summary["rows"] if row.get("vanishes_at_true_root")]
                first = vanishing[0] if vanishing else None
                support_families = []
                if first is not None:
                    support_families = sorted({item["family"] for item in first["support"]})
                sweep_row = {
                    "m": m_value,
                    "t": t_value,
                    "basis_rows": summary["basis_rows"],
                    "basis_cols": summary["basis_cols"],
                    "lll_elapsed": summary["lll_elapsed"],
                    "vanishing_rows": len(vanishing),
                    "first_vanishing_row": None if first is None else first["row"],
                    "first_vanishing_norm_bits": None if first is None else first["norm_bits"],
                    "first_vanishing_support_count": None if first is None else first["support_count"],
                    "first_vanishing_support_families": support_families[:16],
                }
                sweep_rows.append(sweep_row)
                if not args.json:
                    first_text = "-"
                    if first is not None:
                        first_text = (
                            f"row{first['row']} norm={first['norm_bits']} "
                            f"support={first['support_count']}"
                        )
                    print(
                        f"construction={args.construction} m={m_value} t={t_value} "
                        f"x_origin={args.x_origin} x_shift_limit={args.x_shift_limit} "
                        f"dim={summary['basis_rows']}x{summary['basis_cols']} "
                        f"lll={summary['lll_elapsed']} vanish={len(vanishing)} first={first_text}",
                        flush=True,
                    )
        if args.json:
            print(
                json.dumps(
                    {
                        "status": "sweep",
                        "construction": args.construction,
                        "x_shift_limit": args.x_shift_limit,
                        "x_origin": args.x_origin,
                        "prime_bits": args.prime_bits,
                        "low_free_width": args.low_free_width,
                        "low_free_offset": model.low_free_offset,
                        "variables": variables,
                        "results": sweep_rows,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        return 0

    if args.mode == "prune":
        try:
            full_summary, source = first_vanishing_source(
                model,
                args.m,
                args.t,
                root_values,
                args.construction,
                args.x_shift_limit,
                args.x_origin,
            )
        except NoVanishingRow as ex:
            full_summary = ex.summary
            result = {
                "status": "no_source_vanishing_row",
                "construction": args.construction,
                "x_shift_limit": args.x_shift_limit,
                "x_origin": args.x_origin,
                "prime_bits": args.prime_bits,
                "low_free_width": args.low_free_width,
                "low_free_offset": model.low_free_offset,
                "m": args.m,
                "t": args.t,
                "full_basis_rows": full_summary["basis_rows"],
                "full_basis_cols": full_summary["basis_cols"],
            }
            if args.json:
                print(json.dumps(result, indent=2, sort_keys=True))
            else:
                print(
                    f"status=no_source_vanishing_row construction={args.construction} "
                    f"x_origin={args.x_origin} x_shift_limit={args.x_shift_limit} "
                    f"m={args.m} t={args.t} "
                    f"dim={full_summary['basis_rows']}x{full_summary['basis_cols']}"
                )
            return 1

        keep_families = {
            item["family"]
            for item in source["support"]
            if abs(item["coeff"]).bit_length() >= args.min_support_coeff_bits
        }
        if not keep_families:
            raise SystemExit("support coefficient threshold removed every source family")
        pruned_summary = profile_model(
            model,
            args.m,
            args.t,
            root_values,
            0,
            keep_families,
            args.construction,
            args.x_shift_limit,
            args.x_origin,
        )
        pruned_vanishing = [row for row in pruned_summary["rows"] if row.get("vanishes_at_true_root")]
        result = {
            "status": "prune",
            "construction": args.construction,
            "x_shift_limit": args.x_shift_limit,
            "x_origin": args.x_origin,
                "prime_bits": args.prime_bits,
                "low_free_width": args.low_free_width,
                "low_free_offset": model.low_free_offset,
                "m": args.m,
            "t": args.t,
            "variables": variables,
            "source": {
                "basis_rows": full_summary["basis_rows"],
                "basis_cols": full_summary["basis_cols"],
                "lll_elapsed": full_summary["lll_elapsed"],
                "vanishing_rows": len([row for row in full_summary["rows"] if row.get("vanishes_at_true_root")]),
                "first_vanishing_row": source["row"],
                "first_vanishing_norm_bits": source["norm_bits"],
                "support_count": source["support_count"],
                "min_support_coeff_bits": args.min_support_coeff_bits,
            },
            "pruned": {
                "basis_rows": pruned_summary["basis_rows"],
                "basis_cols": pruned_summary["basis_cols"],
                "lll_elapsed": pruned_summary["lll_elapsed"],
                "retained_families": len(keep_families),
                "vanishing_rows": len(pruned_vanishing),
                "first_vanishing_row": None if not pruned_vanishing else pruned_vanishing[0]["row"],
                "first_vanishing_norm_bits": None if not pruned_vanishing else pruned_vanishing[0]["norm_bits"],
            },
            "retained_family_sample": sorted(keep_families)[:24],
        }
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True))
        else:
            print(
                f"status=prune construction={args.construction} m={args.m} t={args.t} "
                f"x_origin={args.x_origin} x_shift_limit={args.x_shift_limit} "
                f"source_dim={full_summary['basis_rows']}x{full_summary['basis_cols']} "
                f"source_vanish={result['source']['vanishing_rows']} source_row={source['row']} "
                f"min_coeff_bits={args.min_support_coeff_bits} retained={len(keep_families)} "
                f"pruned_dim={pruned_summary['basis_rows']}x{pruned_summary['basis_cols']} "
                f"pruned_vanish={len(pruned_vanishing)} "
                f"first_pruned={result['pruned']['first_vanishing_row']}"
            )
        return 0

    if args.mode == "drop-sweep":
        try:
            full_summary, source = first_vanishing_source(
                model,
                args.m,
                args.t,
                root_values,
                args.construction,
                args.x_shift_limit,
                args.x_origin,
            )
        except NoVanishingRow:
            print(
                f"status=no_source_vanishing_row construction={args.construction} "
                f"x_origin={args.x_origin} x_shift_limit={args.x_shift_limit} m={args.m} t={args.t}"
            )
            return 1
        source_families = {item["family"] for item in source["support"]}
        axes = [axis.strip() for axis in args.drop_axes.split(",") if axis.strip()]
        tests = []
        for axis in axes:
            if axis == "power":
                values = sorted({item["power"] for item in source["support"]})
                for value in values:
                    drop = {item["family"] for item in source["support"] if item["power"] == value}
                    tests.append((axis, str(value), drop))
            elif axis in ("nscale", "n_scale_power"):
                values = sorted({item["n_scale_power"] for item in source["support"]})
                for value in values:
                    drop = {item["family"] for item in source["support"] if item["n_scale_power"] == value}
                    tests.append(("nscale", str(value), drop))
            elif axis in ("y", "y_shift"):
                values = sorted({item["y_shift"] for item in source["support"]})
                for value in values:
                    drop = {item["family"] for item in source["support"] if item["y_shift"] == value}
                    tests.append(("y_shift", str(value), drop))
            elif axis in ("x", "x_shift"):
                values = sorted({item["x_shift"] for item in source["support"]})
                for value in values:
                    drop = {item["family"] for item in source["support"] if item["x_shift"] == value}
                    tests.append(("x_shift", str(value), drop))
            else:
                raise SystemExit(f"unknown drop axis: {axis}")

        rows = []
        for axis, value, drop in tests:
            keep_families = source_families - drop
            if not keep_families:
                rows.append(
                    {
                        "axis": axis,
                        "value": value,
                        "dropped": len(drop),
                        "retained": 0,
                        "vanishing_rows": 0,
                        "first_vanishing_row": None,
                        "status": "empty",
                    }
                )
                continue
            pruned_summary = profile_model(
                model,
                args.m,
                args.t,
                root_values,
                0,
                keep_families,
                args.construction,
                args.x_shift_limit,
                args.x_origin,
            )
            pruned_vanishing = [row for row in pruned_summary["rows"] if row.get("vanishes_at_true_root")]
            row = {
                "axis": axis,
                "value": value,
                "dropped": len(drop),
                "retained": len(keep_families),
                "basis_rows": pruned_summary["basis_rows"],
                "basis_cols": pruned_summary["basis_cols"],
                "lll_elapsed": pruned_summary["lll_elapsed"],
                "vanishing_rows": len(pruned_vanishing),
                "first_vanishing_row": None if not pruned_vanishing else pruned_vanishing[0]["row"],
                "status": "retained" if pruned_vanishing else "lost",
            }
            rows.append(row)
            if not args.json:
                print(
                    f"axis={axis} value={value} drop={len(drop)} retained={len(keep_families)} "
                    f"dim={row['basis_rows']}x{row['basis_cols']} vanish={len(pruned_vanishing)} "
                    f"status={row['status']}",
                    flush=True,
                )
        if args.json:
            print(
                json.dumps(
                    {
                        "status": "drop-sweep",
                        "construction": args.construction,
                        "x_shift_limit": args.x_shift_limit,
                        "x_origin": args.x_origin,
                        "prime_bits": args.prime_bits,
                        "low_free_width": args.low_free_width,
                        "low_free_offset": model.low_free_offset,
                        "m": args.m,
                        "t": args.t,
                        "source": {
                            "basis_rows": full_summary["basis_rows"],
                            "basis_cols": full_summary["basis_cols"],
                            "first_vanishing_row": source["row"],
                            "support_count": source["support_count"],
                        },
                        "results": rows,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        return 0

    summary = profile_model(
        model,
        args.m,
        args.t,
        root_values,
        args.rows,
        construction=args.construction,
        x_shift_limit=args.x_shift_limit,
        x_origin=args.x_origin,
    )
    summary.update(
        {
            "prime_bits": args.prime_bits,
            "low_free_width": args.low_free_width,
            "low_free_offset": model.low_free_offset,
            "variables": variables,
        }
    )
    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(
            "status=profile "
            f"construction={args.construction} "
            f"x_origin={args.x_origin} x_shift_limit={args.x_shift_limit} "
            f"prime_bits={args.prime_bits} low_free_width={args.low_free_width} "
            f"low_free_offset={model.low_free_offset} "
            f"m={args.m} t={args.t} dim={summary['basis_rows']}x{summary['basis_cols']} "
            f"lll_elapsed={summary['lll_elapsed']:.3f}"
        )
        for row in summary["rows"]:
            if row.get("status") == "not_divisible":
                print(f"row={row['row']} status=not_divisible")
                continue
            supports = ",".join(
                f"{item['power']}:{item['x_shift']}:{item['y_shift']}:{item['coeff']}"
                for item in row["support"][:6]
            )
            print(
                f"row={row['row']} norm_bits={row['norm_bits']} terms={row['terms']} "
                f"degree={row['degree']} vanishes={row['vanishes_at_true_root']} "
                f"value_bits={row['value_bits_at_true_root']} support_count={row['support_count']} "
                f"support={supports}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
