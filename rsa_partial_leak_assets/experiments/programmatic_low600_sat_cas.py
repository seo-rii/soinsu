#!/usr/bin/env python3
"""PySAT outer loop for partial-low600 SAT+CAS experiments.

This is the first scaffold for `ct07_programmatic_low600_sat_cas`.
It can either build the existing full multiplication CNF for one 4-bit edge
candidate or build a smaller prefix CNF for `p*q == N mod 2^k`.  The prefix
mode is the closer Ajani-Bright-style path: get a low-bit SAT projection, call
the partial-low600 cuso oracle, then optionally learn a blocking clause.

By default cuso `soft_no_root` results are not learned as hard clauses.  Use
`--learn-soft-no-root` only for controlled experiments after a soundness gate.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

ASSET_ROOT = Path(__file__).resolve().parents[1]
for path in (
    ASSET_ROOT / ".py-site",
    ASSET_ROOT / ".sage-site",
    ASSET_ROOT / "src",
    ASSET_ROOT / "experiments",
):
    if path.exists() and str(path) not in sys.path:
        sys.path.insert(0, str(path))

from solve7_main import LEAK, MASK, N, testp  # noqa: E402
from sat_mul_cnf import SatBuilderCNF, build as build_sat  # noqa: E402
from low600_partial_cuso import (  # noqa: E402
    SHAPES,
    Segment,
    build_model,
    model_record,
    print_record,
    run_cuso,
    run_leaf_certify,
)
from challenge_focus_template_hm import (  # noqa: E402
    expand_model_by_small_variables,
    run_one_model as run_focus_template_model,
)


def model_values(raw_model: list[int], variable_count: int) -> list[bool]:
    values = [False] * (variable_count + 1)
    for lit in raw_model:
        if lit > 0:
            values[lit] = True
    return values


def bit_from_sat(bit: int, pvar_by_bit: dict[int, int], values: list[bool], base: int) -> int:
    var = pvar_by_bit.get(bit)
    if var is None:
        return (base >> bit) & 1
    return 1 if values[var] else 0


def fixed_id_from_sat(
    shape: str,
    pvar_by_bit: dict[int, int],
    values: list[bool],
    base: int,
) -> tuple[int, list[dict[str, int]]]:
    fixed_id = 0
    shift = 0
    fixed = []
    for seg in SHAPES[shape]:
        value = 0
        for offset in range(seg.width):
            value |= bit_from_sat(seg.start + offset, pvar_by_bit, values, base) << offset
        fixed_id |= value << shift
        shift += seg.width
        fixed.append({"start": seg.start, "width": seg.width, "value": value})
    return fixed_id, fixed


def blocking_clause_for_projection(
    shape: str,
    pvar_by_bit: dict[int, int],
    values: list[bool],
    base: int,
) -> list[int]:
    clause = []
    for seg in SHAPES[shape]:
        for offset in range(seg.width):
            bit = seg.start + offset
            var = pvar_by_bit.get(bit)
            if var is None:
                continue
            value = bit_from_sat(bit, pvar_by_bit, values, base)
            clause.append(-var if value else var)
    return clause


def blocking_clause_for_fixed_segments(
    fixed: list[dict[str, int]],
    pvar_by_bit: dict[int, int],
    base: int,
) -> list[int]:
    clause = []
    for item in fixed:
        start = int(item["start"])
        width = int(item["width"])
        value = int(item["value"])
        for offset in range(width):
            bit = start + offset
            var = pvar_by_bit.get(bit)
            bit_value = (value >> offset) & 1
            if var is None:
                if bit_value != ((base >> bit) & 1):
                    return []
                continue
            clause.append(-var if bit_value else var)
    return clause


def projection_phase_literals(
    shape: str,
    pvar_by_bit: dict[int, int],
    policy: str,
    seed: int,
    iteration: int,
) -> list[int]:
    if policy == "default":
        return []
    rng = random.Random(seed + iteration)
    phases = []
    ordinal = 0
    for seg in SHAPES[shape]:
        for offset in range(seg.width):
            bit = seg.start + offset
            var = pvar_by_bit.get(bit)
            if var is None:
                continue
            if policy == "zero":
                value = False
            elif policy == "one":
                value = True
            elif policy == "alternate":
                value = bool(ordinal & 1)
            elif policy == "random":
                value = bool(rng.getrandbits(1))
            else:
                raise ValueError(f"unknown projection phase policy {policy!r}")
            phases.append(var if value else -var)
            ordinal += 1
    return phases


def candidate_p(meta, values: list[bool]) -> int:
    p = meta["base"]
    for bit, var in meta["pvars"]:
        if values[var]:
            p |= 1 << bit
        else:
            p &= ~(1 << bit)
    return p


def print_json(record: dict, as_json: bool) -> None:
    if as_json:
        print(json.dumps(record, sort_keys=True), flush=True)
    else:
        print(" ".join(f"{key}={value}" for key, value in record.items()), flush=True)


def append_jsonl(path_text: str | None, record: dict) -> None:
    if not path_text:
        return
    path = Path(path_text)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as out:
        out.write(json.dumps(record, sort_keys=True) + "\n")


def projection_key(shape: str, fixed_id: int) -> str:
    return f"{shape}:{fixed_id}"


def focus_template_args(args):
    return argparse.Namespace(
        mode="run",
        m=args.focus_m,
        t=args.focus_t,
        construction="x-limited",
        x_shift_limit=args.focus_x_shift_limit,
        template_families=args.focus_template_families,
        rows=args.focus_rows,
        delta=args.focus_delta,
        precision=args.focus_precision,
        recover=args.focus_recover,
        prime_start=args.focus_prime_start,
        prime_count=args.focus_prime_count,
        max_candidates=args.focus_max_candidates,
        stop_on_empty_prime=args.focus_stop_on_empty_prime,
        verbose_recover=args.focus_verbose_recover,
    )


def run_focus_template(model, args, queue_index: int | None = None) -> list[dict]:
    origins = ("lower", "upper") if args.focus_x_origin == "both" else (args.focus_x_origin,)
    focus_args = focus_template_args(args)
    records = []
    for sub_model, brute_assignments in expand_model_by_small_variables(
        model,
        args.focus_brute_small_vars,
    ):
        for origin in origins:
            record = run_focus_template_model(sub_model, focus_args, origin, queue_index)
            if brute_assignments:
                record["brute_assignments"] = brute_assignments
            records.append(record)
    return records


def projection_queue_record(
    args,
    iteration: int,
    fixed_id: int,
    fixed: list[dict[str, int]],
    blocking_clause_len: int,
) -> dict:
    return {
        "record_type": "prefix_projection",
        "version": 1,
        "key": projection_key(args.shape, fixed_id),
        "shape": args.shape,
        "fixed_id": fixed_id,
        "fixed": fixed,
        "fixed_bits": sum(item["width"] for item in fixed),
        "blocking_clause_len": blocking_clause_len,
        "prefix_bits": args.prefix_bits,
        "solver": args.solver,
        "iteration": iteration,
        "projection_phase": args.projection_phase,
        "phase_seed": args.phase_seed,
    }


def iter_queue_records(path_text: str):
    path = Path(path_text)
    with path.open("r", encoding="utf-8") as src:
        for lineno, line in enumerate(src, 1):
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            if record.get("record_type") != "prefix_projection":
                continue
            yield lineno, record


def apply_learned_clauses(builder, args, pvar_by_bit: dict[int, int], base: int) -> None:
    loaded = 0
    added = 0
    literals = 0
    for raw_path in args.load_learned_jsonl:
        path = Path(raw_path)
        with path.open("r", encoding="utf-8") as src:
            for line in src:
                line = line.strip()
                if not line:
                    continue
                record = json.loads(line)
                if record.get("record_type") != "learned_projection_clause":
                    continue
                loaded += 1
                clause = blocking_clause_for_fixed_segments(record.get("fixed", []), pvar_by_bit, base)
                if not clause:
                    continue
                builder.s.add_clause(clause)
                added += 1
                literals += len(clause)
    if args.load_learned_jsonl:
        print_json(
            {
                "event": "loaded_learned_clauses",
                "files": len(args.load_learned_jsonl),
                "records": loaded,
                "clauses": added,
                "literals": literals,
            },
            args.json,
        )


def build_prefix_sat(prefix_bits: int, solver_name: str):
    if not 1 <= prefix_bits <= 1024:
        raise ValueError("prefix_bits must be in 1..1024")
    builder = SatBuilderCNF(solver_name)
    pbits = []
    qbits = []
    pvars = []
    for bit in range(prefix_bits):
        if (MASK >> bit) & 1:
            pbits.append(bool((LEAK >> bit) & 1))
        else:
            var = builder.new()
            pbits.append(var)
            pvars.append((bit, var))
    for bit in range(prefix_bits):
        var = builder.new()
        qbits.append(var)
        if bit == 0:
            builder.eqconst(var, True)

    cols = [[] for _ in range(prefix_bits + 1)]
    const = [0] * (prefix_bits + 1)
    for i, pbit in enumerate(pbits):
        if pbit is False:
            continue
        for j, qbit in enumerate(qbits[: prefix_bits - i]):
            col = i + j
            if pbit is True:
                term = qbit
            elif qbit is True:
                term = pbit
            elif qbit is False:
                continue
            else:
                term = builder.and2(pbit, qbit)
            if term is True:
                const[col] += 1
            elif term is not False:
                cols[col].append(term)

    for bit in range(prefix_bits):
        if const[bit] >= 2:
            const[bit + 1] += const[bit] // 2
            const[bit] %= 2
        column = cols[bit]
        if const[bit] == 1:
            column.append(True)
            const[bit] = 0
        while len(column) > 2:
            a = column.pop()
            b = column.pop()
            c = column.pop()
            sum_bit, carry = builder.full_adder(a, b, c)
            if sum_bit is True:
                column.append(True)
            elif sum_bit is not False:
                column.append(sum_bit)
            if carry is True:
                const[bit + 1] += 1
            elif carry is not False:
                cols[bit + 1].append(carry)
            true_count = sum(1 for item in column if item is True)
            if true_count >= 2:
                column[:] = [item for item in column if item is not True]
                const[bit + 1] += true_count // 2
                if true_count & 1:
                    column.append(True)
        cols[bit] = column

    carry = False
    for bit in range(prefix_bits):
        if const[bit] >= 2:
            const[bit + 1] += const[bit] // 2
            const[bit] %= 2
        column = cols[bit]
        if const[bit] == 1:
            column = column + [True]
        while len(column) > 2:
            a = column.pop()
            b = column.pop()
            c = column.pop()
            sum_bit, carry_out = builder.full_adder(a, b, c)
            if sum_bit is not False:
                column.append(sum_bit)
            if carry_out is True:
                const[bit + 1] += 1
            elif carry_out is not False:
                cols[bit + 1].append(carry_out)
        a = column[0] if len(column) > 0 else False
        b = column[1] if len(column) > 1 else False
        sum_bit, carry = builder.full_adder(a, b, carry)
        builder.eqconst(sum_bit, (N >> bit) & 1)

    return builder, {"base": LEAK, "pbits": pbits, "qbits": qbits, "pvars": pvars}


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["self-test", "loop", "prefix-loop", "replay-queue"], default="loop")
    ap.add_argument("--shape", choices=sorted(SHAPES), default="B")
    ap.add_argument("--edge-low", type=int, default=0)
    ap.add_argument("--edge-high", type=int, default=0)
    ap.add_argument("--prefix-bits", type=int, default=600)
    ap.add_argument("--solver", default="cadical153")
    ap.add_argument("--iterations", type=int, default=1)
    ap.add_argument("--conf-budget", type=int, default=100000)
    ap.add_argument("--oracle", choices=["dry-run", "cuso", "focus-template", "leaf-certify"], default="dry-run")
    ap.add_argument("--queue-jsonl", help="append prefix projection records or replay this queue")
    ap.add_argument("--results-jsonl", help="append replay/oracle result records")
    ap.add_argument("--load-learned-jsonl", action="append", default=[])
    ap.add_argument("--queue-start", type=int, default=0, help="0-based replay record start")
    ap.add_argument("--queue-limit", type=int, help="maximum replay records")
    ap.add_argument(
        "--projection-phase",
        choices=["default", "zero", "one", "alternate", "random"],
        default="default",
    )
    ap.add_argument("--phase-seed", type=int, default=20260706)
    ap.add_argument("--learn-soft-no-root", action="store_true")
    ap.add_argument("--upper-variable", action="append", default=[])
    ap.add_argument("--upper-low-variables", action="store_true")
    ap.add_argument("--upper-all-variables", action="store_true")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--cuso-log", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--cuso-graph", choices=["auto", "on", "off"], default="auto")
    ap.add_argument("--cuso-no-intermediate", action="store_true")
    ap.add_argument("--cuso-allow-partial", action="store_true")
    ap.add_argument("--leaf-max-completions", type=int, default=16)
    ap.add_argument("--leaf-start", type=int, default=0)
    ap.add_argument("--leaf-stop", type=int)
    ap.add_argument("--focus-x-origin", choices=["lower", "upper", "both"], default="both")
    ap.add_argument("--focus-m", type=int, default=11)
    ap.add_argument("--focus-t", type=int, default=3)
    ap.add_argument("--focus-x-shift-limit", type=int, default=1)
    ap.add_argument("--focus-template-families", default="")
    ap.add_argument("--focus-rows", type=int, default=13)
    ap.add_argument("--focus-delta", type=float, default=0.99)
    ap.add_argument("--focus-precision", type=int, default=192)
    ap.add_argument("--focus-recover", choices=["none", "resultant", "groebner", "both"], default="both")
    ap.add_argument("--focus-prime-start", type=int, default=257)
    ap.add_argument("--focus-prime-count", type=int, default=56)
    ap.add_argument("--focus-max-candidates", type=int, default=2000)
    ap.add_argument(
        "--focus-brute-small-vars",
        type=int,
        default=0,
        metavar="WIDTH",
        help="brute-force focus-template variables up to WIDTH bits until two variables remain",
    )
    ap.add_argument("--focus-no-stop-on-empty-prime", dest="focus_stop_on_empty_prime", action="store_false")
    ap.set_defaults(focus_stop_on_empty_prime=True)
    ap.add_argument("--focus-verbose-recover", action="store_true")
    return ap.parse_args(argv)


def self_test() -> None:
    pvar_by_bit = {bit: bit - 264 for bit in range(265, 349)}
    values = [False] * 90
    for bit, var in pvar_by_bit.items():
        values[var] = bool(bit & 1)
    fixed_id, fixed = fixed_id_from_sat("C", pvar_by_bit, values, LEAK)
    expected = sum(((bit & 1) << (bit - 265)) for bit in range(265, 349))
    assert fixed_id == expected, (fixed_id, expected)
    assert fixed == [{"start": 265, "width": 84, "value": expected}]
    clause = blocking_clause_for_projection("C", pvar_by_bit, values, LEAK)
    assert len(clause) == 84
    fixed_clause = blocking_clause_for_fixed_segments(
        [{"start": 265, "width": 2, "value": 1}],
        {265: 10, 266: 11},
        LEAK,
    )
    assert fixed_clause == [-10, 11]

    pvar_by_bit = {bit: bit - 361 for bit in range(362, 420)}
    values = [False] * 64
    for bit, var in pvar_by_bit.items():
        values[var] = bool((bit >> 1) & 1)
    clause = blocking_clause_for_projection("B", pvar_by_bit, values, LEAK)
    assert len(clause) == 58
    phases = projection_phase_literals("B", pvar_by_bit, "zero", 1, 0)
    assert len(phases) == 58 and all(lit < 0 for lit in phases)
    phases = projection_phase_literals("B", pvar_by_bit, "one", 1, 0)
    assert len(phases) == 58 and all(lit > 0 for lit in phases)
    queue_args = argparse.Namespace(
        shape="B",
        prefix_bits=600,
        solver="cadical153",
        projection_phase="random",
        phase_seed=3,
    )
    queued = projection_queue_record(
        queue_args,
        0,
        7,
        [{"start": 150, "width": 4, "value": 7}, {"start": 362, "width": 58, "value": 0}],
        62,
    )
    assert queued["key"] == "B:7"
    assert queued["fixed_bits"] == 62
    builder, meta = build_prefix_sat(16, "cadical153")
    assert builder.nv > 0
    assert builder.s.solve()
    values = model_values(builder.s.get_model(), builder.nv)
    p_low = 0
    q_low = 0
    for bit, item in enumerate(meta["pbits"]):
        value = item if isinstance(item, bool) else values[item]
        if value:
            p_low |= 1 << bit
    for bit, item in enumerate(meta["qbits"]):
        value = item if isinstance(item, bool) else values[item]
        if value:
            q_low |= 1 << bit
    assert (p_low * q_low - N) % (1 << 16) == 0
    focus_args = argparse.Namespace(
        focus_m=11,
        focus_t=3,
        focus_x_shift_limit=1,
        focus_template_families="",
        focus_rows=13,
        focus_delta=0.99,
        focus_precision=192,
        focus_recover="none",
        focus_prime_start=257,
        focus_prime_count=2,
        focus_max_candidates=4,
        focus_brute_small_vars=0,
        focus_stop_on_empty_prime=True,
        focus_verbose_recover=False,
        focus_x_origin="both",
    )
    focus_model = build_model("E_150_4_265_84_362_32", 0, [])
    focus_records = run_focus_template(focus_model, focus_args)
    assert [record["x_origin"] for record in focus_records] == ["lower", "upper"]
    assert all(record["basis_rows"] == 13 and record["basis_cols"] == 15 for record in focus_records)
    print("self-test ok")


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.mode == "self-test":
        self_test()
        return 0

    if args.mode == "replay-queue":
        if not args.queue_jsonl:
            print("--queue-jsonl is required for replay-queue", file=sys.stderr)
            return 2
        processed = 0
        exit_status = 1
        for index, (lineno, queued) in enumerate(iter_queue_records(args.queue_jsonl)):
            if index < args.queue_start:
                continue
            if args.queue_limit is not None and processed >= args.queue_limit:
                break
            processed += 1
            shape = queued["shape"]
            fixed_id = int(queued["fixed_id"])
            model = build_model(
                shape,
                fixed_id,
                [],
                args.upper_variable,
                args.upper_low_variables,
                args.upper_all_variables,
            )
            replay_record = {
                "event": "queue_replay",
                "queue_line": lineno,
                "queue_index": index,
                "key": queued.get("key", projection_key(shape, fixed_id)),
                "shape": shape,
                "fixed_id": fixed_id,
                "oracle": args.oracle,
            }
            print_json(replay_record, args.json)
            append_jsonl(args.results_jsonl, replay_record)
            if args.oracle == "dry-run":
                record = model_record(model, "dry_run")
                print_record(record, args.json)
                append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **replay_record, **record})
                exit_status = 0
                continue
            if args.oracle == "focus-template":
                for record in run_focus_template(model, args, index):
                    print_json(record, args.json)
                    append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **replay_record, **record})
                    if record["status"] == "factor":
                        return 0
                exit_status = 1
                continue
            if args.oracle == "leaf-certify":
                leaf_records, summary = run_leaf_certify(model, args)
                for record in leaf_records:
                    print_json(record, args.json)
                    append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **replay_record, **record})
                    if record["status"] == "factor":
                        return 0
                print_json(summary, args.json)
                append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **replay_record, **summary})
                if summary["status"] == "factor":
                    return 0
                exit_status = 1
                continue
            status, roots_count, oracle_elapsed = run_cuso(model, args)
            record = model_record(model, status, roots_count, oracle_elapsed)
            print_record(record, args.json)
            append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **replay_record, **record})
            if status == "factor":
                return 0
            if status in ("candidate", "soft_no_root"):
                exit_status = 1
        if processed == 0:
            print_json({"event": "queue_empty", "queue": args.queue_jsonl}, args.json)
            return 1
        return exit_status

    if args.mode == "prefix-loop":
        required_prefix = max((seg.stop for seg in SHAPES[args.shape]), default=0)
        if args.prefix_bits < required_prefix:
            print_json(
                {
                    "event": "prefix_too_small",
                    "shape": args.shape,
                    "prefix_bits": args.prefix_bits,
                    "required_prefix_bits": required_prefix,
                },
                args.json,
            )
            return 2
        t0 = time.time()
        builder, meta = build_prefix_sat(args.prefix_bits, args.solver)
        pvar_by_bit = {bit: var for bit, var in meta["pvars"]}
        apply_learned_clauses(builder, args, pvar_by_bit, meta["base"])
        print_json(
            {
                "event": "built_prefix_sat",
                "solver": args.solver,
                "prefix_bits": args.prefix_bits,
                "vars": builder.nv,
                "clauses": builder.clauses,
                "elapsed": round(time.time() - t0, 3),
            },
            args.json,
        )
        for iteration in range(args.iterations):
            phases = projection_phase_literals(
                args.shape,
                pvar_by_bit,
                args.projection_phase,
                args.phase_seed,
                iteration,
            )
            if phases:
                builder.s.set_phases(phases)
                print_json(
                    {
                        "event": "set_projection_phases",
                        "iteration": iteration,
                        "policy": args.projection_phase,
                        "lits": len(phases),
                    },
                    args.json,
                )
            builder.s.conf_budget(args.conf_budget)
            st = time.time()
            sat = builder.s.solve_limited(expect_interrupt=True)
            elapsed = time.time() - st
            print_json(
                {
                    "event": "prefix_sat_solve",
                    "iteration": iteration,
                    "status": sat,
                    "elapsed": round(elapsed, 3),
                },
                args.json,
            )
            if sat is not True:
                return 1 if sat is None else 0
            values = model_values(builder.s.get_model(), builder.nv)
            fixed_id, fixed = fixed_id_from_sat(args.shape, pvar_by_bit, values, meta["base"])
            clause = blocking_clause_for_projection(args.shape, pvar_by_bit, values, meta["base"])
            print_json(
                {
                    "event": "prefix_projection",
                    "iteration": iteration,
                    "shape": args.shape,
                    "fixed_id": fixed_id,
                    "fixed": fixed,
                    "blocking_clause_len": len(clause),
                },
                args.json,
            )
            queued = projection_queue_record(
                args,
                iteration,
                fixed_id,
                fixed,
                len(clause),
            )
            append_jsonl(args.queue_jsonl, queued)
            if args.queue_jsonl:
                print_json(
                    {
                        "event": "queued_projection",
                        "iteration": iteration,
                        "key": queued["key"],
                        "queue": args.queue_jsonl,
                    },
                    args.json,
                )
            oracle_model = build_model(
                args.shape,
                fixed_id,
                [],
                args.upper_variable,
                args.upper_low_variables,
                args.upper_all_variables,
            )
            if args.oracle == "dry-run":
                print_record(model_record(oracle_model, "dry_run"), args.json)
                if clause:
                    builder.s.add_clause(clause)
                    print_json(
                        {
                            "event": "blocked_projection_for_sampling",
                            "iteration": iteration,
                            "lits": len(clause),
                        },
                        args.json,
                    )
                continue
            if args.oracle == "focus-template":
                records = run_focus_template(oracle_model, args)
                for record in records:
                    print_json(record, args.json)
                    append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **record})
                    if record["status"] == "factor":
                        return 0
                print_json(
                    {
                        "event": "stop_without_hard_clause",
                        "iteration": iteration,
                        "oracle_status": ",".join(record["status"] for record in records),
                    },
                    args.json,
                )
                return 1
            if args.oracle == "leaf-certify":
                leaf_records, summary = run_leaf_certify(oracle_model, args)
                for record in leaf_records:
                    print_json(record, args.json)
                    append_jsonl(
                        args.results_jsonl,
                        {
                            "record_type": "oracle_result",
                            "event": "prefix_leaf_result",
                            "iteration": iteration,
                            "oracle": args.oracle,
                            **record,
                        },
                    )
                    if record["status"] == "factor":
                        return 0
                print_json(summary, args.json)
                append_jsonl(
                    args.results_jsonl,
                    {
                        "record_type": "oracle_result",
                        "event": "prefix_leaf_summary",
                        "iteration": iteration,
                        "oracle": args.oracle,
                        **summary,
                    },
                )
                if summary["status"] == "factor":
                    return 0
                if summary["status"] == "leaf_exhausted_soft_no_root" and args.learn_soft_no_root:
                    if not clause:
                        event = {"event": "skip_empty_clause", "iteration": iteration}
                        print_json(event, args.json)
                        append_jsonl(args.results_jsonl, event)
                        return 1
                    builder.s.add_clause(clause)
                    event = {
                        "event": "learned_clause",
                        "iteration": iteration,
                        "kind": "leaf_exhausted_soft_no_root_promoted",
                        "lits": len(clause),
                        "shape": args.shape,
                        "fixed_id": fixed_id,
                    }
                    print_json(event, args.json)
                    append_jsonl(args.results_jsonl, event)
                    continue
                event = {
                    "event": "stop_without_hard_clause",
                    "iteration": iteration,
                    "oracle_status": summary["status"],
                    "shape": args.shape,
                    "fixed_id": fixed_id,
                }
                print_json(event, args.json)
                append_jsonl(args.results_jsonl, event)
                return 1
            status, roots_count, oracle_elapsed = run_cuso(oracle_model, args)
            record = model_record(oracle_model, status, roots_count, oracle_elapsed)
            print_record(record, args.json)
            append_jsonl(
                args.results_jsonl,
                {
                    "record_type": "oracle_result",
                    "event": "prefix_oracle_result",
                    "iteration": iteration,
                    "oracle": args.oracle,
                    **record,
                },
            )
            if status == "factor":
                return 0
            if status == "soft_no_root" and args.learn_soft_no_root:
                if not clause:
                    event = {"event": "skip_empty_clause", "iteration": iteration}
                    print_json(event, args.json)
                    append_jsonl(args.results_jsonl, event)
                    return 1
                builder.s.add_clause(clause)
                event = {
                    "event": "learned_clause",
                    "iteration": iteration,
                    "kind": "soft_no_root_promoted",
                    "lits": len(clause),
                    "shape": args.shape,
                    "fixed_id": fixed_id,
                }
                print_json(event, args.json)
                append_jsonl(args.results_jsonl, event)
                continue
            event = {
                "event": "stop_without_hard_clause",
                "iteration": iteration,
                "oracle_status": status,
                "shape": args.shape,
                "fixed_id": fixed_id,
            }
            print_json(event, args.json)
            append_jsonl(args.results_jsonl, event)
            return 1
        return 0

    t0 = time.time()
    builder, meta = build_sat(N, MASK, LEAK, args.edge_low, args.edge_high, args.solver)
    pvar_by_bit = {bit: var for bit, var in meta["pvars"]}
    apply_learned_clauses(builder, args, pvar_by_bit, meta["base"])
    print_json(
        {
            "event": "built_sat",
            "solver": args.solver,
            "edge_low": args.edge_low,
            "edge_high": args.edge_high,
            "vars": builder.nv,
            "clauses": builder.clauses,
            "q_prefix_bits": meta["pref"],
            "elapsed": round(time.time() - t0, 3),
        },
        args.json,
    )

    for iteration in range(args.iterations):
        phases = projection_phase_literals(
            args.shape,
            pvar_by_bit,
            args.projection_phase,
            args.phase_seed,
            iteration,
        )
        if phases:
            builder.s.set_phases(phases)
            print_json(
                {
                    "event": "set_projection_phases",
                    "iteration": iteration,
                    "policy": args.projection_phase,
                    "lits": len(phases),
                },
                args.json,
            )
        builder.s.conf_budget(args.conf_budget)
        st = time.time()
        sat = builder.s.solve_limited(expect_interrupt=True)
        elapsed = time.time() - st
        print_json(
            {
                "event": "sat_solve",
                "iteration": iteration,
                "status": sat,
                "elapsed": round(elapsed, 3),
            },
            args.json,
        )
        if sat is not True:
            return 1 if sat is None else 0

        values = model_values(builder.s.get_model(), builder.nv)
        p = candidate_p(meta, values)
        if testp(p):
            return 0

        fixed_id, fixed = fixed_id_from_sat(args.shape, pvar_by_bit, values, meta["base"])
        clause = blocking_clause_for_projection(args.shape, pvar_by_bit, values, meta["base"])
        record = {
            "event": "projection",
            "iteration": iteration,
            "shape": args.shape,
            "fixed_id": fixed_id,
            "fixed": fixed,
            "blocking_clause_len": len(clause),
        }
        print_json(record, args.json)

        oracle_model = build_model(
            args.shape,
            fixed_id,
            [],
            args.upper_variable,
            args.upper_low_variables,
            args.upper_all_variables,
        )
        if args.oracle == "dry-run":
            print_record(model_record(oracle_model, "dry_run"), args.json)
            return 0

        if args.oracle == "focus-template":
            records = run_focus_template(oracle_model, args)
            for record in records:
                print_json(record, args.json)
                append_jsonl(args.results_jsonl, {"record_type": "oracle_result", **record})
                if record["status"] == "factor":
                    return 0
            print_json(
                {
                    "event": "stop_without_hard_clause",
                    "iteration": iteration,
                    "oracle_status": ",".join(record["status"] for record in records),
                },
                args.json,
            )
            return 1

        if args.oracle == "leaf-certify":
            leaf_records, summary = run_leaf_certify(oracle_model, args)
            for record in leaf_records:
                print_json(record, args.json)
                append_jsonl(
                    args.results_jsonl,
                    {
                        "record_type": "oracle_result",
                        "event": "leaf_result",
                        "iteration": iteration,
                        "oracle": args.oracle,
                        **record,
                    },
                )
                if record["status"] == "factor":
                    return 0
            print_json(summary, args.json)
            append_jsonl(
                args.results_jsonl,
                {
                    "record_type": "oracle_result",
                    "event": "leaf_summary",
                    "iteration": iteration,
                    "oracle": args.oracle,
                    **summary,
                },
            )
            if summary["status"] == "factor":
                return 0
            if summary["status"] == "leaf_exhausted_soft_no_root" and args.learn_soft_no_root:
                if not clause:
                    event = {"event": "skip_empty_clause", "iteration": iteration}
                    print_json(event, args.json)
                    append_jsonl(args.results_jsonl, event)
                    return 1
                builder.s.add_clause(clause)
                event = {
                    "event": "learned_clause",
                    "iteration": iteration,
                    "kind": "leaf_exhausted_soft_no_root_promoted",
                    "lits": len(clause),
                    "shape": args.shape,
                    "fixed_id": fixed_id,
                }
                print_json(event, args.json)
                append_jsonl(args.results_jsonl, event)
                continue
            event = {
                "event": "stop_without_hard_clause",
                "iteration": iteration,
                "oracle_status": summary["status"],
                "shape": args.shape,
                "fixed_id": fixed_id,
            }
            print_json(event, args.json)
            append_jsonl(args.results_jsonl, event)
            return 1

        status, roots_count, oracle_elapsed = run_cuso(oracle_model, args)
        record = model_record(oracle_model, status, roots_count, oracle_elapsed)
        print_record(record, args.json)
        append_jsonl(
            args.results_jsonl,
            {
                "record_type": "oracle_result",
                "event": "oracle_result",
                "iteration": iteration,
                "oracle": args.oracle,
                **record,
            },
        )
        if status == "factor":
            return 0
        if status == "soft_no_root" and args.learn_soft_no_root:
            if not clause:
                event = {"event": "skip_empty_clause", "iteration": iteration}
                print_json(event, args.json)
                append_jsonl(args.results_jsonl, event)
                return 1
            builder.s.add_clause(clause)
            event = {
                "event": "learned_clause",
                "iteration": iteration,
                "kind": "soft_no_root_promoted",
                "lits": len(clause),
                "shape": args.shape,
                "fixed_id": fixed_id,
            }
            print_json(event, args.json)
            append_jsonl(args.results_jsonl, event)
            continue
        event = {
            "event": "stop_without_hard_clause",
            "iteration": iteration,
            "oracle_status": status,
            "shape": args.shape,
            "fixed_id": fixed_id,
        }
        print_json(event, args.json)
        append_jsonl(args.results_jsonl, event)
        return 1
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
