#!/usr/bin/env python3
"""Generate a compact challenge CNF for Ajani--Bright's programmatic MapleSAT.

The first 2,048 DIMACS variables are deliberately reserved for every bit of
``p`` and ``q`` (least-significant bit first).  Known bits are still emitted as
unit clauses, while the multiplication circuit substitutes their Boolean
values to keep the circuit substantially smaller.  This gives the MapleSAT
Coppersmith callback the contiguous bit-variable layout it expects without
paying for gates involving constants.

The generated DIMACS header intentionally contains placeholder counts.  The
MapleSAT parser grows its variable table while reading clauses and only warns
about a header mismatch.  Streaming this way avoids retaining millions of
clauses in memory or writing a second uncompressed temporary file.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
import sys
import time
from pathlib import Path

ASSET_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ASSET_ROOT / "src"))

from solve7_main import LEAK, MASK, N  # noqa: E402


P_BITS = 1024
Q_BITS = 1024
P_LSB_VAR = 1
P_MSB_VAR = P_BITS
Q_LSB_VAR = P_BITS + 1
Q_MSB_VAR = P_BITS + Q_BITS


class StreamingCNF:
    def __init__(self, output: Path):
        self.output = output
        if output.suffix == ".gz":
            self._stream = gzip.open(output, "wt", encoding="ascii", compresslevel=3)
        else:
            self._stream = output.open("w", encoding="ascii")
        self.nv = Q_MSB_VAR
        self.clauses = 0
        self._stream.write(f"c {N}\n")
        self._stream.write(f"c {P_LSB_VAR}\n")
        self._stream.write(f"c {P_MSB_VAR}\n")
        self._stream.write(f"c {Q_LSB_VAR}\n")
        self._stream.write(f"c {Q_MSB_VAR}\n")
        self._stream.write("c p and q variables are contiguous and LSB-first\n")
        self._stream.write("p cnf 1 0\n")

    def close(self) -> None:
        self._stream.close()

    def new(self) -> int:
        self.nv += 1
        return self.nv

    def clause(self, literals) -> None:
        out = []
        seen = set()
        for literal in literals:
            if literal is True:
                return
            if literal is False:
                continue
            literal = int(literal)
            if -literal in seen:
                return
            if literal not in seen:
                seen.add(literal)
                out.append(literal)
        if not out:
            raise ValueError("attempted to emit an empty clause")
        self._stream.write(" ".join(map(str, out)))
        self._stream.write(" 0\n")
        self.clauses += 1

    def unit(self, literal, value: bool = True) -> None:
        if literal is True:
            if not value:
                raise ValueError("inconsistent true constant")
            return
        if literal is False:
            if value:
                raise ValueError("inconsistent false constant")
            return
        self.clause([literal if value else -literal])

    def xor_bit(self, arguments, constant: bool = False):
        variables = []
        parity = constant
        for argument in arguments:
            if argument is False:
                continue
            if argument is True:
                parity = not parity
            else:
                variables.append(int(argument))
        if not variables:
            return parity
        if len(variables) == 1:
            return -variables[0] if parity else variables[0]
        output = self.new()
        all_variables = [output, *variables]
        for bits in itertools.product((False, True), repeat=len(all_variables)):
            if (sum(bits) & 1) != parity:
                self.clause(
                    -variable if bit else variable
                    for variable, bit in zip(all_variables, bits)
                )
        return output

    def and2(self, left, right):
        if left is False or right is False:
            return False
        if left is True:
            return right
        if right is True:
            return left
        if left == right:
            return left
        if left == -right:
            return False
        output = self.new()
        self.clause([-output, left])
        self.clause([-output, right])
        self.clause([output, -left, -right])
        return output

    def or2(self, left, right):
        if left is True or right is True:
            return True
        if left is False:
            return right
        if right is False:
            return left
        if left == right:
            return left
        if left == -right:
            return True
        output = self.new()
        self.clause([output, -left])
        self.clause([output, -right])
        self.clause([-output, left, right])
        return output

    def majority3(self, left, middle, right):
        arguments = (left, middle, right)
        true_count = sum(argument is True for argument in arguments)
        variables = [
            argument
            for argument in arguments
            if argument is not True and argument is not False
        ]
        needed = 2 - true_count
        if needed <= 0:
            return True
        if needed > len(variables):
            return False
        if needed == 1:
            output = variables[0]
            for variable in variables[1:]:
                output = self.or2(output, variable)
            return output
        if needed == len(variables):
            output = variables[0]
            for variable in variables[1:]:
                output = self.and2(output, variable)
            return output
        left, middle, right = variables
        output = self.new()
        self.clause([-left, -middle, output])
        self.clause([-left, -right, output])
        self.clause([-middle, -right, output])
        self.clause([left, middle, -output])
        self.clause([left, right, -output])
        self.clause([middle, right, -output])
        return output

    def full_adder(self, left, middle, right):
        return self.xor_bit([left, middle, right]), self.majority3(left, middle, right)

    def eqconst(self, value, bit: int) -> None:
        self.unit(value, bool(bit))


def contiguous_known_lsb_bits(mask: int, total_bits: int) -> int:
    width = 0
    while width < total_bits and ((mask >> width) & 1):
        width += 1
    return width


def common_prefix_interval(lower: int, upper: int, total_bits: int) -> tuple[int, int]:
    differing = lower ^ upper
    width = total_bits - differing.bit_length() if differing else total_bits
    return width, lower >> (total_bits - width) if width else 0


def add_bits(cnf: StreamingCNF, left, right, carry=False):
    """Return the little-endian ripple-carry sum of two bit vectors."""
    width = max(len(left), len(right))
    sums = []
    for bit in range(width):
        left_bit = left[bit] if bit < len(left) else False
        right_bit = right[bit] if bit < len(right) else False
        sum_bit, carry = cnf.full_adder(left_bit, right_bit, carry)
        sums.append(sum_bit)
    sums.append(carry)
    return sums


def wallace_multiply(cnf: StreamingCNF, left, right):
    """Multiply little-endian bit vectors with a Wallace reduction tree."""
    rows = []
    for right_index, right_bit in enumerate(right):
        row = [False] * right_index
        row.extend(cnf.and2(left_bit, right_bit) for left_bit in left)
        rows.append(row)
    if not rows:
        return []
    while len(rows) > 2:
        reduced = []
        for start in range(0, len(rows) - 2, 3):
            group = rows[start : start + 3]
            width = max(map(len, group))
            carries = [False]
            sums = []
            for bit in range(width):
                arguments = [row[bit] if bit < len(row) else False for row in group]
                sum_bit, carry_bit = cnf.full_adder(*arguments)
                sums.append(sum_bit)
                carries.append(carry_bit)
            reduced.extend((carries, sums))
        reduced.extend(rows[(len(rows) // 3) * 3 :])
        rows = reduced
    return rows[0] if len(rows) == 1 else add_bits(cnf, rows[0], rows[1])


def recursive_multiply(cnf: StreamingCNF, left, right, bound: int = 20):
    """Purdom--Sabry recursive (Karatsuba) multiplier."""
    size = max(len(left), len(right))
    if size <= bound:
        return wallace_multiply(cnf, left, right)

    half = (size + 1) // 2
    left_size = len(left)
    right_size = len(right)
    padded_left = list(left) + [False] * (size - left_size)
    padded_right = list(right) + [False] * (size - right_size)
    left_low, left_high = padded_left[:half], padded_left[half:]
    right_low, right_high = padded_right[:half], padded_right[half:]

    if right_size <= half:
        low_product = recursive_multiply(cnf, left_low, right_low, bound)
        high_product = recursive_multiply(cnf, left_high, right_low, bound)
        return add_bits(cnf, [False] * half + high_product, low_product)

    high_product = recursive_multiply(cnf, left_high, right_high, bound)
    low_product = recursive_multiply(cnf, left_low, right_low, bound)
    left_sum = add_bits(cnf, left_high, left_low)
    right_sum = add_bits(cnf, right_high, right_low)
    cross_with_ends = recursive_multiply(cnf, left_sum, right_sum, bound)

    width = max(len(cross_with_ends), len(high_product))
    complemented_high = []
    for bit in range(width):
        value = high_product[bit] if bit < len(high_product) else False
        complemented_high.append(False if value is True else True if value is False else -value)
    cross_plus_low = add_bits(cnf, cross_with_ends, complemented_high, True)
    cnf.unit(cross_plus_low[-1], True)
    cross_plus_low = cross_plus_low[:-1]

    width = max(len(cross_plus_low), len(low_product))
    complemented_low = []
    for bit in range(width):
        value = low_product[bit] if bit < len(low_product) else False
        complemented_low.append(False if value is True else True if value is False else -value)
    cross = add_bits(cnf, cross_plus_low, complemented_low, True)
    cnf.unit(cross[-1], True)
    cross = cross[:-1]

    high_and_cross = add_bits(
        cnf,
        [False] * (2 * half) + high_product,
        [False] * half + cross,
    )
    return add_bits(cnf, high_and_cross, low_product)


def build_known_bits(edge_low: int | None, edge_high: int | None):
    fixed_mask = MASK
    fixed_value = LEAK
    if edge_low is not None:
        fixed_mask |= 0xF << 150
        fixed_value |= edge_low << 150
    if edge_high is not None:
        fixed_mask |= 0xF << 920
        fixed_value |= edge_high << 920

    p_lsb_width = contiguous_known_lsb_bits(fixed_mask, P_BITS)
    modulus = 1 << p_lsb_width
    p_lsb = fixed_value & (modulus - 1)
    q_lsb = (N * pow(p_lsb, -1, modulus)) % modulus

    all_bits = (1 << P_BITS) - 1
    p_min = fixed_value
    p_max = fixed_value | (all_bits & ~fixed_mask)
    q_min = N // p_max
    q_max = N // p_min
    q_prefix_width, q_prefix = common_prefix_interval(q_min, q_max, Q_BITS)

    pbits = []
    qbits = []
    for bit in range(P_BITS):
        variable = P_LSB_VAR + bit
        if (fixed_mask >> bit) & 1:
            pbits.append(bool((fixed_value >> bit) & 1))
        else:
            pbits.append(variable)
    for bit in range(Q_BITS):
        variable = Q_LSB_VAR + bit
        if bit < p_lsb_width:
            qbits.append(bool((q_lsb >> bit) & 1))
        elif bit >= Q_BITS - q_prefix_width:
            qbits.append(bool((q_prefix >> (bit - (Q_BITS - q_prefix_width))) & 1))
        else:
            qbits.append(variable)

    metadata = {
        "edge_low": edge_low,
        "edge_high": edge_high,
        "p_fixed_bits": fixed_mask.bit_count(),
        "p_lsb_width": p_lsb_width,
        "q_lsb_width": p_lsb_width,
        "q_prefix_width": q_prefix_width,
        "p_unknown_below_569": sum(
            1 for bit in range(569) if not ((fixed_mask >> bit) & 1)
        ),
        "p_unknown_below_618": sum(
            1 for bit in range(618) if not ((fixed_mask >> bit) & 1)
        ),
    }
    return fixed_mask, fixed_value, pbits, qbits, metadata


def generate(
    output: Path,
    edge_low: int | None,
    edge_high: int | None,
    multiplier: str = "csa",
) -> dict:
    fixed_mask, fixed_value, pbits, qbits, metadata = build_known_bits(edge_low, edge_high)
    output.parent.mkdir(parents=True, exist_ok=True)
    cnf = StreamingCNF(output)
    started = time.time()
    try:
        for bit in range(P_BITS):
            if (fixed_mask >> bit) & 1:
                cnf.unit(P_LSB_VAR + bit, bool((fixed_value >> bit) & 1))
        for bit, value in enumerate(qbits):
            if value is True or value is False:
                cnf.unit(Q_LSB_VAR + bit, value)

        edge_implications = 0
        if edge_low is None:
            edge_variables = [P_LSB_VAR + bit for bit in range(150, 154)]
            for candidate in range(16):
                guards = [
                    -variable if (candidate >> offset) & 1 else variable
                    for offset, variable in enumerate(edge_variables)
                ]
                _, _, _, candidate_qbits, _ = build_known_bits(candidate, edge_high)
                for bit, value in enumerate(candidate_qbits):
                    if qbits[bit] is True or qbits[bit] is False:
                        continue
                    if value is True or value is False:
                        cnf.clause([*guards, (Q_LSB_VAR + bit) if value else -(Q_LSB_VAR + bit)])
                        edge_implications += 1
        if edge_high is None:
            edge_variables = [P_LSB_VAR + bit for bit in range(920, 924)]
            for candidate in range(16):
                guards = [
                    -variable if (candidate >> offset) & 1 else variable
                    for offset, variable in enumerate(edge_variables)
                ]
                _, _, _, candidate_qbits, _ = build_known_bits(edge_low, candidate)
                for bit, value in enumerate(candidate_qbits):
                    if qbits[bit] is True or qbits[bit] is False:
                        continue
                    if value is True or value is False:
                        cnf.clause([*guards, (Q_LSB_VAR + bit) if value else -(Q_LSB_VAR + bit)])
                        edge_implications += 1

        if multiplier == "recursive":
            product = recursive_multiply(cnf, pbits, qbits)
            for bit, product_bit in enumerate(product):
                cnf.eqconst(product_bit, (N >> bit) & 1 if bit < 2048 else 0)
            terms = 0
        else:
            max_columns = 2060
            columns = [[] for _ in range(max_columns)]
            constants = [0] * max_columns
            terms = 0
            for i, pbit in enumerate(pbits):
                if pbit is False:
                    continue
                for j, qbit in enumerate(qbits):
                    column = i + j
                    if pbit is True:
                        term = qbit
                    elif qbit is True:
                        term = pbit
                    elif qbit is False:
                        continue
                    else:
                        term = cnf.and2(pbit, qbit)
                    if term is True:
                        constants[column] += 1
                    elif term is not False:
                        columns[column].append(term)
                    terms += 1

            for bit in range(max_columns - 1):
                if constants[bit] >= 2:
                    constants[bit + 1] += constants[bit] // 2
                    constants[bit] %= 2
                column = columns[bit]
                if constants[bit] == 1:
                    column.append(True)
                    constants[bit] = 0
                while len(column) > 2:
                    left = column.pop()
                    middle = column.pop()
                    right = column.pop()
                    sum_bit, carry = cnf.full_adder(left, middle, right)
                    if sum_bit is True:
                        column.append(True)
                    elif sum_bit is not False:
                        column.append(sum_bit)
                    if carry is True:
                        constants[bit + 1] += 1
                    elif carry is not False:
                        columns[bit + 1].append(carry)
                    true_count = sum(item is True for item in column)
                    if true_count >= 2:
                        column[:] = [item for item in column if item is not True]
                        constants[bit + 1] += true_count // 2
                        if true_count & 1:
                            column.append(True)
                columns[bit] = column

            carry = False
            for bit in range(max_columns - 1):
                if constants[bit] >= 2:
                    constants[bit + 1] += constants[bit] // 2
                    constants[bit] %= 2
                column = columns[bit]
                if constants[bit] == 1:
                    column = column + [True]
                while len(column) > 2:
                    left = column.pop()
                    middle = column.pop()
                    right = column.pop()
                    sum_bit, carry_out = cnf.full_adder(left, middle, right)
                    if sum_bit is not False:
                        column.append(sum_bit)
                    if carry_out is True:
                        constants[bit + 1] += 1
                    elif carry_out is not False:
                        columns[bit + 1].append(carry_out)
                left = column[0] if column else False
                right = column[1] if len(column) > 1 else False
                sum_bit, carry = cnf.full_adder(left, right, carry)
                cnf.eqconst(sum_bit, (N >> bit) & 1 if bit < 2048 else 0)
            cnf.eqconst(carry, 0)
    finally:
        cnf.close()

    metadata.update(
        {
            "output": str(output),
            "variables": cnf.nv,
            "clauses": cnf.clauses,
            "multiplication_terms": terms,
            "multiplier": multiplier,
            "edge_implication_clauses": edge_implications,
            "elapsed": time.time() - started,
            "compressed_bytes": output.stat().st_size,
        }
    )
    return metadata


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--edge-low", type=int, choices=range(16))
    parser.add_argument("--edge-high", type=int, choices=range(16))
    parser.add_argument("--multiplier", choices=("csa", "recursive"), default="csa")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    record = generate(args.output, args.edge_low, args.edge_high, args.multiplier)
    print(json.dumps(record, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
