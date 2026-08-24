#!/usr/bin/env python3
"""Row-pruned HM runner for challenge partial-low600 projections.

This applies the planted focus-group 13-family template to real challenge
partial-low600 two-variable shapes such as:

    E_150_4_265_84_362_16: x = p[378..419], y = p[600..1023]
    E_150_4_265_84_362_32: x = p[394..419], y = p[600..1023]

It is an experimental fallback lattice runner.  A failure to recover roots is
not a sound no-root certificate.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ASSET_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("DOT_SAGE", str(ASSET_ROOT / ".sage-tmp"))
for path in (ASSET_ROOT / "src", ASSET_ROOT / "experiments"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from fpylll import LLL  # noqa: E402
from low600_partial_cuso import SHAPES, build_model  # noqa: E402
from planted_focus_group_hm import build_basis, row_norm_bits, row_to_poly  # noqa: E402
from root_groebner import recover as recover_groebner, roots_mod_groebner  # noqa: E402
from root_recover import recover_crt, roots_mod_pair  # noqa: E402
from solve7_main import N, testp  # noqa: E402


DEFAULT_TEMPLATE_FAMILIES = {
    "0:0:3",
    "0:1:3",
    "1:0:2",
    "1:1:0:2",
    "1:1:2",
    "2:0:1",
    "2:1:0:1",
    "2:1:1",
    "2:1:1:1",
    "2:2:1",
    "3:0:0",
    "3:1:0",
    "3:1:0:0",
}


def print_record(record: dict, as_json: bool) -> None:
    if as_json:
        print(json.dumps(record, sort_keys=True), flush=True)
        return
    print(" ".join(f"{key}={value}" for key, value in record.items()), flush=True)


def iter_queue(path_text: str):
    with Path(path_text).open("r", encoding="utf-8") as src:
        for index, line in enumerate(src):
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get("record_type") == "prefix_projection":
                yield index, record


def first_empty_prime(polys, primes, recover: str, verbose: bool) -> dict | None:
    checks = 0
    enabled_paths = [
        (path_name, root_fn)
        for path_name, root_fn in (
            ("resultant", roots_mod_pair),
            ("groebner", roots_mod_groebner),
        )
        if recover in (path_name, "both")
    ]
    for prime in primes:
        nonempty_paths = []
        for path_name, root_fn in enabled_paths:
            checks += 1
            pairs = root_fn(polys, prime)
            if verbose:
                print(f"{path_name} prime {prime} pairs {len(pairs)}", flush=True)
            if pairs:
                nonempty_paths.append(path_name)
        if not nonempty_paths:
            return {
                "empty_prime": prime,
                "empty_path": "+".join(path_name for path_name, _ in enabled_paths),
                "prime_checks": checks,
            }
    return None


def expand_model_by_small_variables(model: dict, max_width: int):
    expanded = [(model, [])]
    if max_width <= 0:
        return expanded
    while any(len(item_model["variables"]) > 2 for item_model, _ in expanded):
        next_expanded = []
        changed = False
        for item_model, assignments in expanded:
            variables = item_model["variables"]
            if len(variables) <= 2:
                next_expanded.append((item_model, assignments))
                continue
            brute_seg = next((seg for seg in variables if seg.width <= max_width), None)
            if brute_seg is None:
                next_expanded.append((item_model, assignments))
                continue
            changed = True
            for value in range(1 << brute_seg.width):
                sub_model = dict(item_model)
                sub_model["base"] = item_model["base"] + (value << brute_seg.start)
                sub_model["fixed"] = list(item_model["fixed"]) + [(brute_seg, value)]
                sub_model["variables"] = [seg for seg in variables if seg != brute_seg]
                sub_model["fixed_bits"] = item_model["fixed_bits"] + brute_seg.width
                sub_model["variable_bits"] = item_model["variable_bits"] - brute_seg.width
                sub_model["variable_origins"] = {
                    seg: origin
                    for seg, origin in item_model.get("variable_origins", {}).items()
                    if seg != brute_seg
                }
                next_expanded.append(
                    (
                        sub_model,
                        assignments
                        + [
                            {
                                "start": brute_seg.start,
                                "width": brute_seg.width,
                                "value": value,
                            }
                        ],
                    )
                )
        expanded = next_expanded
        if not changed:
            break
    return expanded


def run_one_model(model: dict, args, x_origin: str, queue_index: int | None = None) -> dict:
    variables = model["variables"]
    if len(variables) != 2:
        raise ValueError(
            "challenge focus-template runner expects exactly two variables; "
            f"{model['shape']} has {len(variables)}.  Use --brute-small-vars "
            "to enumerate small variables first."
        )
    template = {
        family.strip()
        for family in args.template_families.split(",")
        if family.strip()
    } or set(DEFAULT_TEMPLATE_FAMILIES)
    challenge_model = SimpleNamespace(
        base=model["base"],
        n=N,
        variables=variables,
    )
    matrix, monomials, scales, shifts = build_basis(
        challenge_model,
        args.m,
        args.t,
        template,
        args.construction,
        args.x_shift_limit,
        "upper" if x_origin == "upper" else "zero",
    )
    if args.mode == "dry-run":
        return {
            "status": "dry_run",
            "shape": model["shape"],
            "fixed_id": model["fixed_id"],
            "queue_index": queue_index,
            "x_origin": x_origin,
            "m": args.m,
            "t": args.t,
            "template_families": len(template),
            "basis_rows": matrix.nrows,
            "basis_cols": matrix.ncols,
            "variables": [
                {"start": seg.start, "width": seg.width}
                for seg in variables
            ],
        }

    start_time = time.time()
    LLL.reduction(
        matrix,
        delta=args.delta,
        method="proved",
        float_type="mpfr",
        precision=args.precision,
    )
    lll_elapsed = time.time() - start_time

    rows = []
    for row_index in range(matrix.nrows):
        poly = row_to_poly(matrix, row_index, monomials, scales)
        if poly and len(poly) > 1:
            rows.append((row_norm_bits(matrix, row_index), row_index, poly))
    rows.sort(key=lambda item: (item[0], item[1]))
    selected = [poly for _, _, poly in rows[: args.rows]]
    row_stats = [
        {"row": row_index, "norm_bits": round(norm_bits, 3), "terms": len(poly)}
        for norm_bits, row_index, poly in rows[: min(args.rows, 8)]
    ]

    candidates = []
    recover_elapsed = 0.0
    empty_prime_info = None
    if args.recover != "none" and len(selected) >= 2:
        import sympy as sp

        primes = []
        prime = args.prime_start - 1
        for _ in range(args.prime_count):
            prime = int(sp.nextprime(prime))
            primes.append(prime)
        recover_start = time.time()
        if args.stop_on_empty_prime:
            empty_prime_info = first_empty_prime(
                selected,
                primes,
                args.recover,
                args.verbose_recover,
            )
        if empty_prime_info is None:
            if args.recover in ("resultant", "both"):
                candidates.extend(
                    recover_crt(
                        selected,
                        1 << variables[0].width,
                        1 << variables[1].width,
                        primes,
                        verbose=args.verbose_recover,
                    )
                )
            if not candidates and args.recover in ("groebner", "both"):
                candidates.extend(
                    recover_groebner(
                        selected,
                        1 << variables[0].width,
                        1 << variables[1].width,
                        primes,
                        verbose=args.verbose_recover,
                    )
                )
        recover_elapsed = time.time() - recover_start

    tested = 0
    for x_root, y_root in candidates[: args.max_candidates]:
        tested += 1
        x_value = ((1 << variables[0].width) - 1) - x_root if x_origin == "upper" else x_root
        p = model["base"] + (x_value << variables[0].start) + (y_root << variables[1].start)
        if testp(p):
            return {
                "status": "factor",
                "shape": model["shape"],
                "fixed_id": model["fixed_id"],
                "queue_index": queue_index,
                "x_origin": x_origin,
                "m": args.m,
                "t": args.t,
                "basis_rows": matrix.nrows,
                "basis_cols": matrix.ncols,
                "lll_elapsed": round(lll_elapsed, 3),
                "recover_elapsed": round(recover_elapsed, 3),
                "rows": len(selected),
                "candidate_count": len(candidates),
                "tested": tested,
                "p": str(p),
            }

    return {
        "status": "not_found",
        "shape": model["shape"],
        "fixed_id": model["fixed_id"],
        "queue_index": queue_index,
        "x_origin": x_origin,
        "m": args.m,
        "t": args.t,
        "template_families": len(template),
        "basis_rows": matrix.nrows,
        "basis_cols": matrix.ncols,
        "lll_elapsed": round(lll_elapsed, 3),
        "recover_elapsed": round(recover_elapsed, 3),
        "rows": len(selected),
        "row_stats": row_stats,
        "candidate_count": len(candidates),
        "tested": tested,
        "recover": args.recover,
        **({} if empty_prime_info is None else empty_prime_info),
    }


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["self-test", "dry-run", "run"], default="dry-run")
    ap.add_argument("--shape", choices=sorted(SHAPES), default="E_150_4_265_84_362_16")
    ap.add_argument("--fixed-id", type=int, default=0)
    ap.add_argument("--queue-jsonl")
    ap.add_argument("--queue-start", type=int, default=0)
    ap.add_argument("--queue-limit", type=int)
    ap.add_argument("--x-origin", choices=["lower", "upper", "both"], default="both")
    ap.add_argument(
        "--brute-small-vars",
        type=int,
        default=0,
        metavar="WIDTH",
        help="brute-force variables up to WIDTH bits until a two-variable model remains",
    )
    ap.add_argument("--m", type=int, default=11)
    ap.add_argument("--t", type=int, default=3)
    ap.add_argument("--construction", choices=["x-limited"], default="x-limited")
    ap.add_argument("--x-shift-limit", type=int, default=1)
    ap.add_argument("--template-families", default="")
    ap.add_argument("--rows", type=int, default=13)
    ap.add_argument("--delta", type=float, default=0.99)
    ap.add_argument("--precision", type=int, default=192)
    ap.add_argument("--recover", choices=["none", "resultant", "groebner", "both"], default="resultant")
    ap.add_argument("--prime-start", type=int, default=257)
    ap.add_argument("--prime-count", type=int, default=56)
    ap.add_argument("--max-candidates", type=int, default=2000)
    ap.add_argument("--no-stop-on-empty-prime", dest="stop_on_empty_prime", action="store_false")
    ap.set_defaults(stop_on_empty_prime=True)
    ap.add_argument("--verbose-recover", action="store_true")
    ap.add_argument("--json", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.mode == "self-test":
        assert first_empty_prime([{(1, 0): 1}, {(0, 1): 1}], [257], "both", False) is None
        empty = first_empty_prime([{(0, 0): 1}, {(1, 0): 1}], [257], "both", False)
        assert empty is not None and empty["empty_prime"] == 257, empty
        for shape in ("E_150_4_265_84_362_16", "E_150_4_265_84_362_32"):
            model = build_model(shape, 0, [])
            for origin in ("lower", "upper"):
                record = run_one_model(model, argparse.Namespace(**{**vars(args), "mode": "dry-run"}), origin)
                assert record["basis_rows"] == 13, record
                assert record["basis_cols"] == 15, record
        shape_c = build_model("C", 0, [])
        expanded_c = expand_model_by_small_variables(shape_c, 4)
        assert len(expanded_c) == 16, len(expanded_c)
        assert all(len(item_model["variables"]) == 2 for item_model, _ in expanded_c)
        print("self-test ok")
        return 0

    origins = ("lower", "upper") if args.x_origin == "both" else (args.x_origin,)
    exit_status = 1
    records = []
    if args.queue_jsonl:
        processed = 0
        for queue_index, queued in iter_queue(args.queue_jsonl):
            if queue_index < args.queue_start:
                continue
            if args.queue_limit is not None and processed >= args.queue_limit:
                break
            processed += 1
            model = build_model(queued["shape"], int(queued["fixed_id"]), [])
            for sub_model, brute_assignments in expand_model_by_small_variables(
                model,
                args.brute_small_vars,
            ):
                for origin in origins:
                    record = run_one_model(sub_model, args, origin, queue_index)
                    if brute_assignments:
                        record["brute_assignments"] = brute_assignments
                    records.append(record)
                    print_record(record, args.json)
                    if record["status"] == "factor":
                        return 0
                    if record["status"] == "dry_run":
                        exit_status = 0
        if processed == 0:
            print_record({"status": "queue_empty", "queue_jsonl": args.queue_jsonl}, args.json)
            return 1
    else:
        model = build_model(args.shape, args.fixed_id, [])
        for sub_model, brute_assignments in expand_model_by_small_variables(
            model,
            args.brute_small_vars,
        ):
            for origin in origins:
                record = run_one_model(sub_model, args, origin)
                if brute_assignments:
                    record["brute_assignments"] = brute_assignments
                records.append(record)
                print_record(record, args.json)
                if record["status"] == "factor":
                    return 0
                if record["status"] == "dry_run":
                    exit_status = 0
    if any(record["status"] == "not_found" for record in records):
        return 1
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
