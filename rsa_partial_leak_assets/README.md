# RSA partial-prime-bit leak challenge #7 assets

This package contains the reproducible work assets for the 2026 cryptanalysis contest problem #7.
The original problem statement is in `problem/problem_7_partial_factorization.pdf`.

## Current status

The final RSA factorization/plaintext has **not** been recovered in this package yet.  The package preserves the analysis, scripts, and recommended external run path so that the remaining heavy lattice/SAT jobs can be continued on a Sage/cuso/flatter-capable machine.

## Main idea

The leaked mask fixes 672 of the 1024 bits of the secret prime `p`; 352 bits remain unknown in seven blocks:

```text
bits 150..153    len 4
bits 265..348    len 84
bits 362..419    len 58
bits 600..668    len 69
bits 682..768    len 87
bits 784..829    len 46
bits 920..923    len 4
```

The preferred attack brute-forces the two 4-bit blocks and groups the remaining unknowns into two variables:

```text
x = p[265..419]      155-bit grouped variable
y = p[600..829]      230-bit grouped variable
p = p0 + 2^265*x + 2^600*y
```

Then each of the 256 `(low, high)` candidates is solved as a bivariate Coppersmith/Herrmann-May problem modulo an unknown divisor `p` of the known multiple `N`.

## Fast path: Sage + cuso

Install the dependencies described in `env/setup_external.md`, then run:

```bash
sage -python src/solve7_main.py --mode analyze
sage -python src/solve7_main.py --mode cuso --a 0 --b 256 | tee logs/cuso_full.log
```

To split across workers, candidate id is `cid = 16*high + low`:

```bash
bash scripts/run_cuso_range.sh 0 64
bash scripts/run_cuso_range.sh 64 128
bash scripts/run_cuso_range.sh 128 192
bash scripts/run_cuso_range.sh 192 256
```

## Local fallback

The fallback is slower and was not strong enough in the current container, but is useful for sanity checks:

```bash
python3 -m pip install -r env/requirements-local.txt
python3 src/solve7_main.py --mode local --m 8 --t 3 --lead y --a 0 --b 1
```

## Decryption once a factor is found

```bash
python3 src/decrypt_with_factor.py <p-as-hex-or-decimal>
# or
sage -python src/solve7_main.py --mode decrypt --phex <p-as-hex>
```

Successful output prints `FOUND`, `p`, `q`, `m.hex`, and `m` as a Python bytes literal.

## Directory map

```text
problem/        original PDF problem statement
src/            compact solver, inspection, and decryption helpers
scripts/        shell runners for external machines
env/            setup notes, requirements, and a fallback Dockerfile
docs/           detailed work report in DOCX and Markdown
experiments/    intermediate scripts/logs from the working session
logs/           runtime log destination
```

## References

- Herrmann and May, "Solving Linear Equations Modulo Divisors: On Factoring Given Any Bits", ASIACRYPT 2008.
- Ajani and Bright, "SAT and Lattice Reduction for Integer Factorization", ISSAC 2024.
- `cuso`: "Solving Multivariate Coppersmith Problems" EUROCRYPT 2025 artifact.
- kionactf/coppersmith for alternative Herrmann-May style tooling.
