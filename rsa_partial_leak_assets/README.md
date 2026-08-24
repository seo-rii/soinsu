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

Two cuso paths should be tried before spending more time on SAT-ledger
enumeration.

The legacy grouped path brute-forces the two 4-bit blocks and groups the
remaining unknowns into two variables:

```text
x = p[265..419]      155-bit grouped variable
y = p[600..829]      230-bit grouped variable
p = p0 + 2^265*x + 2^600*y
```

Then each of the 256 `(low, high)` candidates is solved as a bivariate
Coppersmith/Herrmann-May problem modulo an unknown divisor `p` of the known
multiple `N`.

The split-block path keeps the exact unknown block structure.  With
`--cuso-split-brute-small-edges`, the two 4-bit edge blocks are brute-forced and
the remaining five blocks become independent cuso variables:

```text
edge brute force: 150:4, 920:4
cuso variables:   265:84, 362:58, 600:69, 682:87, 784:46
```

Both main cuso modes also support origin variants.  `--cuso-upper-low-variables`
represents low-side variables as `u' = 2^width - 1 - u` while leaving variables
starting at bit 600 or higher in lower-origin form.  Use this lower/upper-low
pair before treating a grouped or split no-root result as shape evidence.

The next experimental line is `ct07_partial_low600_cuso_broad_clause`.  The
goal is to stop treating low600 p-Coppersmith as a fully fixed leaf check and
instead use it as a broader pruning oracle: fix only part of
`150:4,265:84,362:58`, leave the remaining low holes as variables, and group
all high bits as `z600 = p[600..1023]`.

First smoke these shapes:

```text
Shape B:
  fixed:    150:4, 362:58
  variable: 265:84, z600:424
  clause size: 62
  variable mass: 508 bits

Shape C:
  fixed:    265:84
  variable: 150:4, 362:58, z600:424
  clause size: 84
  variable mass: 486 bits

Shape D_265_64_362_8:
  fixed:    265:64, 362:8
  variable: 150:4, 329:20, 370:50, z600:424
  clause size: 72
  variable mass: 498 bits

Shape D_265_48_362_16:
  fixed:    265:48, 362:16
  variable: 150:4, 313:36, 378:42, z600:424
  clause size: 64
  variable mass: 506 bits

Shape E_150_4_265_84_362_16:
  fixed:    150:4, 265:84, 362:16
  variable: 378:42, z600:424
  clause size: 104
  variable mass: 466 bits

Shape E_150_4_265_84_362_32:
  fixed:    150:4, 265:84, 362:32
  variable: 394:26, z600:424
  clause size: 120
  variable mass: 450 bits

Shape F_150_4_265_84_362_50:
  fixed:    150:4, 265:84, 362:50
  variable: 412:8, z600:424
  clause size: 138
  variable mass: 432 bits

Shape F_150_4_265_84_362_54:
  fixed:    150:4, 265:84, 362:54
  variable: 416:4, z600:424
  clause size: 142
  variable mass: 428 bits

Shape F_150_4_265_84_362_58:
  fixed:    150:4, 265:84, 362:58
  variable: z600:424
  clause size: 146
  variable mass: 424 bits

Shape A:
  fixed:    362:58
  variable: 150:4, 265:84, z600:424
  clause size: 58
  variable mass: 512 bits
```

Classify cuso results as `factor`, `candidate`, `soft_no_root`, or
`hard_no_root`.  Do not translate multivariate cuso no-root results into SAT
clauses until a soundness gate is explicit, such as planted true-branch
retention, a stronger oracle check, or another deterministic proof.

The broader literature-backed plan has four named experiment lanes:

```text
ct07_programmatic_low600_sat_cas
  Ajani-Bright-style solver loop:
  SAT model -> partial_low600 cuso -> learned clause for hard_no_root only.
  Scaffold: experiments/programmatic_low600_sat_cas.py.

ct07_cuso_mixed_shape_search
  cuso automatic shift-selection over grouped, exact, mixed, and partial shapes.

ct07_focus_group_hm
  downscaled planted HM/fpylll tests to identify useful lattice rows and prune
  fallback shift bases.

ct07_cocert_clause_minimization
  treat no-root rows as co-certificates and optimize branch coverage per clause.
```

Use these mixed cuso shapes before deciding that grouped or exact is the only
viable model:

```text
S0 grouped 2-var:
  [265:155], [600:230]
  bound mass 385

S1 exact 5-var:
  [265:84], [362:58], [600:69], [682:87], [784:46]
  bound mass 344

S2 low exact + high grouped:
  [265:84], [362:58], [600:230]
  bound mass 372

S3 low grouped + high exact:
  [265:155], [600:69], [682:87], [784:46]
  bound mass 357

S4 low exact + high mixed:
  [265:84], [362:58], [600:169], [784:46]
  bound mass 357

S5 partial low600:
  fixed 58-84 low bits, variables include z600:424
```

## Fast path: Sage + cuso

Install the dependencies described in `env/setup_external.md`, then run:

```bash
mkdir -p .sage-tmp .sage-site
DOT_SAGE="$PWD/.sage-tmp" sage -pip install --target "$PWD/.sage-site" \
  git+https://github.com/keeganryan/cuso.git

export DOT_SAGE="${DOT_SAGE:-$PWD/.sage-tmp}"
export PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}"
sage -python src/solve7_main.py --mode analyze
mkdir -p logs/ct07_cuso

# Smoke the grouped and split cuso modes first.  These smoke commands may
# return "not found" for a single candidate; inspect the logs for cuso import,
# graph, lattice, or root-search errors.
bash scripts/smoke_cuso_modes.sh logs/ct07_cuso/smoke
CUSO_UPPER_LOW_VARIABLES=1 \
  bash scripts/smoke_cuso_modes.sh logs/ct07_cuso/smoke_upper_low

# Full grouped 2-variable sweep.
bash scripts/run_cuso_grouped_8way.sh logs/ct07_cuso/grouped
CUSO_UPPER_LOW_VARIABLES=1 \
  bash scripts/run_cuso_grouped_8way.sh logs/ct07_cuso/grouped_upper_low

# Full split exact-block sweep with 4-bit edge brute force.
bash scripts/run_cuso_split_edges_8way.sh logs/ct07_cuso/split_edges
CUSO_UPPER_LOW_VARIABLES=1 \
  bash scripts/run_cuso_split_edges_8way.sh logs/ct07_cuso/split_edges_upper_low

# Partial-low600 broad-clause smoke.  The first steps run under normal Python;
# the shape cuso steps require Sage+cuso.
bash scripts/smoke_partial_low600_cuso.sh logs/ct07_cuso/partial_low600_smoke
UPPER_LOW_VARIABLES=1 \
  bash scripts/smoke_partial_low600_cuso.sh logs/ct07_cuso/partial_low600_upper_low_smoke

# Mixed grouped/exact/hybrid shape smoke.  This compares S0-S4 for the same
# edge candidate id before full candidate sweeps.
bash scripts/smoke_mixed_shape_cuso.sh logs/ct07_cuso/mixed_shape_smoke
UPPER_LOW_VARIABLES=1 \
  bash scripts/smoke_mixed_shape_cuso.sh logs/ct07_cuso/mixed_shape_upper_low_smoke

# Lower/upper origin matrix for the practical mixed and partial shapes.
bash scripts/run_cuso_origin_matrix.sh logs/ct07_cuso/origin_matrix

# Focused cuso option sweep for the remaining practical shapes.  This compares
# graph on/off, no-intermediate, and allow-partial on S0, S2, and partial B.
bash scripts/run_cuso_option_sweep.sh logs/ct07_cuso/option_sweep

# Downscaled planted true-branch smoke.  This checks that S0/S2/partial-B style
# models can recover a known factor on a small instance before using no-root
# observations as SAT pruning evidence.
bash scripts/smoke_planted_cuso.sh logs/ct07_cuso/planted_smoke

# Configurable planted scale/option sweep.  Use this to narrow the smallest
# scale and cuso option set where a planted true branch is retained.
PLANTED_BITS_LIST=64 PLANTED_SHAPES=S0_grouped_2var \
  PLANTED_OPTIONS='graph_on graph_off_no_intermediate graph_off_allow_partial' \
  CUSO_TIMEOUT=12s \
  bash scripts/run_planted_cuso_scale_sweep.sh logs/ct07_cuso/planted_scale_s0_64

# Planted low-free-width boundary sweep.  This starts at the known-low
# high-tail-only model and opens a small suffix of the low grouped block.
LOW_FREE_WIDTHS=0,1,2,4 CUSO_TIMEOUT=15s \
  bash scripts/run_planted_lowfree_sweep.sh logs/ct07_cuso/planted_lowfree_64

# Focused option sweep for the first failing low-free-width case.
LOW_FREE_WIDTHS=1 CUSO_TIMEOUT=20s \
  bash scripts/run_planted_lowfree_option_sweep.sh logs/ct07_cuso/planted_lowfree1_options

# Focus-group HM profiler for the same planted two-variable low-free model.
python3 experiments/planted_focus_group_hm.py --mode self-test --prime-bits 64
python3 experiments/planted_focus_group_hm.py \
  --mode profile --prime-bits 64 --low-free-width 1 --m 4 --t 2 --rows 12
python3 experiments/planted_focus_group_hm.py \
  --mode sweep --prime-bits 96 --low-free-width 1 --m-values 10-12 --t-values 3-4 \
  --construction x-limited --x-shift-limit 1
python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep --prime-bits 96 --seeds 20260705-20260714 \
  --low-free-width 1 --m-values 11 --t-values 3 \
  --construction x-limited --x-shift-limit 1
python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep --prime-bits 96 --seeds 20260705-20260714 \
  --low-free-width 1 --m-values 11 --t-values 3 \
  --construction x-limited --x-shift-limit 1 --x-origin upper
python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep --prime-bits 96 --seeds 20260705-20260714 \
  --low-free-width 1 --m-values 11 --t-values 3 \
  --construction x-limited --x-shift-limit 1 \
  --template-families '0:0:3,0:1:3,1:0:2,1:1:0:2,1:1:2,2:0:1,2:1:0:1,2:1:1,2:1:1:1,2:2:1,3:0:0,3:1:0,3:1:0:0'
python3 experiments/planted_focus_group_hm.py \
  --mode sweep --prime-bits 64 --low-free-width 1 --m-values 12 --t-values 3-4 \
  --construction total-degree
python3 experiments/planted_focus_group_hm.py \
  --mode prune --prime-bits 64 --low-free-width 1 --m 12 --t 3
python3 experiments/planted_focus_group_hm.py \
  --mode prune --prime-bits 64 --low-free-width 1 --m 12 --t 3 \
  --min-support-coeff-bits 80
python3 experiments/planted_focus_group_hm.py \
  --mode drop-sweep --prime-bits 64 --low-free-width 1 --m 12 --t 3 \
  --drop-axes power,nscale
bash scripts/run_planted_focus_group_sweep.sh logs/ct07_cuso/planted_focus_group_64
PLANTED_BITS_LIST='64 80 96' \
  bash scripts/run_planted_focus_group_scale_sweep.sh logs/ct07_cuso/planted_focus_scale
```

Main grouped lower-origin candidate coverage has also been extended with
bounded per-candidate runs:

```text
latest run: CUSO_TIMEOUT=45s CUSO_GRAPH=off bash scripts/run_cuso_grouped_candidates.sh 40 48
logs:
  logs/ct07_cuso/grouped_candidates_24_32_45s_20260706/
  logs/ct07_cuso/grouped_candidates_32_40_45s_20260706/
  logs/ct07_cuso/grouped_candidates_40_48_45s_20260706/
coverage now recorded:
  lower-origin grouped cid 0..47
  upper-low grouped cid 0..15
cid 24..39 summary:
  timeout 16/16
  multiplicity 3 reached for all 16
  shift polynomials 70 for all 16
  integer relations 0 for all 16
  factor/root not found
cid 40..47 summary:
  timeout 8/8
  multiplicity 3 and shift polynomials 70 for cids 40 and 42..47
  cid 41 timed out at multiplicity 2 with 21 shift polynomials
  integer relations 0 for all 8
  factor/root not found
```

These grouped runs are direct factor-recovery attempts, but the 45-second
budget is only enough to record early cuso progress.  They are not no-root
certificates.

The first SAT-in-loop scaffold is intentionally conservative:

```bash
mkdir -p .py-site
python3 -m pip install --target .py-site 'python-sat[pblib,aiger]>=1.8.dev13'
PYTHONPATH="$PWD/.py-site${PYTHONPATH:+:$PYTHONPATH}" \
  python3 experiments/programmatic_low600_sat_cas.py --mode self-test
PYTHONPATH="$PWD/.py-site${PYTHONPATH:+:$PYTHONPATH}" \
  python3 experiments/programmatic_low600_sat_cas.py \
  --oracle dry-run --shape B --edge-low 0 --edge-high 0 --iterations 1 --json
PYTHONPATH="$PWD/.py-site:$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}" \
  bash scripts/smoke_programmatic_low600_sat_cas.sh \
  logs/ct07_cuso/programmatic_sat_cas_smoke
PROJECTION_PHASE=random PHASE_SEED=7 \
  bash scripts/smoke_programmatic_low600_sat_cas.sh \
  logs/ct07_cuso/programmatic_sat_cas_random_phase_smoke
PROJECTION_PHASE=random PHASE_SEED=11 ITERATIONS=2 \
  bash scripts/smoke_programmatic_low600_sat_cas.sh \
  logs/ct07_cuso/programmatic_sat_cas_queue_smoke
SHAPES='B C' PHASE_SEEDS='1 3 7 11' ITERATIONS=2 \
  bash scripts/run_programmatic_low600_queue_collect.sh \
  logs/ct07_cuso/programmatic_sat_cas_queue_collect
RUN_CUSO=1 CUSO_TIMEOUT=30s CUSO_REPLAY_LIMIT=1 \
  bash scripts/smoke_programmatic_low600_sat_cas.sh \
  logs/ct07_cuso/programmatic_sat_cas_queue_cuso_smoke
CUSO_TIMEOUT=60s QUEUE_LIMIT=8 \
  bash scripts/run_programmatic_low600_queue_cuso.sh \
  logs/ct07_cuso/programmatic_sat_cas_queue_smoke/projections.jsonl \
  logs/ct07_cuso/programmatic_sat_cas_queue_cuso_batch
```

Do not pass `--learn-soft-no-root` for challenge pruning until the selected
partial-low600 shape has a sound no-root gate.  Without that flag, the scaffold
builds the SAT projection and oracle call path but stops before adding an
unsafe learned clause.

On this machine, the first scaffold dry-run for edge `(low=0, high=0)` built
the existing multiplication CNF in 26.598s:

```text
vars: 1,088,050
clauses: 6,820,270
q prefix bits: 192
```

With `--conf-budget 1000`, the SAT call returned `unknown` after 4.991s, so no
SAT projection was available yet.  The next SAT-in-loop step is to tune the
model-generation budget/assumptions before invoking cuso inside the loop.

The newer prefix-SAT path avoids waiting for a full factor model.  It builds
`p*q == N mod 2^600` and samples low600 projections directly:

```text
Shape B prefix 600: build 12.563s, 244,510 vars, 1,550,074 clauses
  projection 0: 0.251s, 62-literal block
  projection 1: 0.091s, 62-literal block

Shape C prefix 600: build 11.899s, 244,510 vars, 1,550,074 clauses
  projection 0: 0.271s, 84-literal block
  projection 1: 0.097s, 84-literal block
```

A 45-second integrated Sage smoke for Shape B prefix 600 reached cuso
graph-off multiplicity 3, rank 75, and 0 integer relations before timeout.
That verifies the intended path now exists:

```text
prefix SAT model -> Shape B/C low600 projection -> partial-low600 cuso oracle
```

This is still not a sound learned-clause loop.  It is the first working
programmatic path where SAT supplies broad low-bit projections without needing
a full factor model.

Prefix-SAT also supports projection phase hints.  This avoids repeatedly
starting from the solver's default all-ones projection:

```text
Shape B phase zero:   fixed_id 0
Shape B phase one:    fixed_id 4611686018427387903
Shape B random seed1: fixed_id 768626895441383966
Shape C random seed7: two distinct 84-literal projections
```

Use `--projection-phase random --phase-seed <n>` or the smoke runner
environment variables `PROJECTION_PHASE=random PHASE_SEED=<n>` when collecting
SAT-supplied low600 projections for oracle experiments.

The smoke runner now writes a reusable projection queue:

```text
logs/.../projections.jsonl
logs/.../results.jsonl
```

In a random-phase queue smoke with `PHASE_SEED=11`, it collected four records:

```text
Shape B: two 62-literal prefix projections
Shape C: two 84-literal prefix projections
Replay: four dry-run oracle models reconstructed from queue
```

For broader collection without cuso, use:

```bash
SHAPES='B C' PHASE_SEEDS='1 3 7 11' ITERATIONS=2 \
  bash scripts/run_programmatic_low600_queue_collect.sh \
  logs/ct07_cuso/programmatic_sat_cas_queue_collect
```

This appends Shape/phase/seed projections to `projections.jsonl` and summarizes
the queue.  A small verification run with `PHASE_SEEDS=17` and `ITERATIONS=1`
collected two unique records:

```text
Shape B: one 62-literal projection
Shape C: one 84-literal projection
duplicates: 0
```

The same queue can be replayed through cuso:

```bash
DOT_SAGE="$PWD/.sage-tmp" timeout 30s \
  sage -python experiments/programmatic_low600_sat_cas.py \
  --mode replay-queue \
  --queue-jsonl logs/ct07_cuso/programmatic_sat_cas_queue_smoke/projections.jsonl \
  --queue-limit 1 \
  --oracle cuso \
  --cuso-graph off \
  --cuso-log INFO \
  --results-jsonl logs/ct07_cuso/programmatic_sat_cas_queue_smoke/cuso_results.jsonl \
  --json
```

The first bounded replay reached cuso graph-off rank 75 with 0 integer
relations before timeout.  This is still evidence of plumbing and cost only,
not a no-root certificate.

The same `replay-queue` path can now run the row-pruned focus-template HM
filter without Sage:

```bash
python3 experiments/programmatic_low600_sat_cas.py \
  --mode replay-queue \
  --queue-jsonl logs/ct07_cuso/programmatic_queue_collect_e_more_20260706/projections.jsonl \
  --queue-limit 2 \
  --oracle focus-template \
  --focus-x-origin both \
  --focus-recover both \
  --results-jsonl logs/ct07_cuso/programmatic_focus_template_replay_smoke_20260706/results.jsonl \
  --json
```

A two-record smoke wrote six JSONL records: two `queue_replay` events and four
focus-template oracle results.  All four oracle results were `not_found` with
`empty_prime=257`; this confirms the programmatic queue can now compare cuso
and focus-template filters through the same replay interface.  These are still
soft heuristic filter results, not hard learned clauses.

For a whole queue, use the focus-template batch runner:

```bash
bash scripts/run_programmatic_low600_queue_focus_template.sh \
  logs/ct07_cuso/programmatic_queue_collect_e_more_20260706/projections.jsonl \
  logs/ct07_cuso/programmatic_focus_template_batch_e_more_20260706
```

The first 16-record E-shape replay produced 32 oracle results because each
projection was tested with lower and upper origins.  All 32 were `not_found`
with `empty_prime=257`, `candidate_count_total=0`, `recover_avg=0.745s`, and
`recover_max=1.475s`.  This makes the focus-template filter cheap enough for
queue-wide soft scoring, but it is still not a hard no-root certificate.

For multiple queue records, prefer the per-record batch runner:

```bash
CUSO_TIMEOUT=60s QUEUE_LIMIT=8 \
  bash scripts/run_programmatic_low600_queue_cuso.sh \
  logs/ct07_cuso/programmatic_sat_cas_queue_smoke/projections.jsonl \
  logs/ct07_cuso/programmatic_sat_cas_queue_cuso_batch
```

It runs one queue index per Sage process and writes one log per record, so a
timeout does not hide which projection reached which cuso stage.  A 20-second
single-record smoke produced:

```text
queue 0, Shape B: timeout, multiplicity 2, shift 24, rank 64, integer relations 0
```

A broader queue run on 2026-07-06 used Shape B/C, random phases, seeds
`1 3 7 11`, and two iterations per seed:

```text
records: 16
unique keys: 16
Shape B: 8 projections, 62-literal fixed block
Shape C: 8 projections, 84-literal fixed block
duplicates: 0
queue: logs/ct07_cuso/programmatic_queue_collect_broad_20260706/projections.jsonl
```

The first four Shape B queue records were replayed with graph-off and a
90-second per-record timeout:

```text
queue 0..3, Shape B lower-origin:
  timeout, multiplicity 3, shift 75, rank 128, integer relations 0
```

Two upper-low-origin replays of the same Shape B records produced the same
summary.  `--cuso-allow-partial` also matched the lower-origin baseline.  With
`--cuso-no-intermediate`, queue 0 jumped to multiplicity 4 and shift 157 within
the same 90-second budget, but still produced 0 integer relations.  A
180-second rerun of the same queue record and option advanced to the
multiplicity-4 lattice at rank 269, again with 0 integer relations before
timeout:

```text
queue 0, Shape B, --cuso-no-intermediate, 180s:
  timeout, multiplicity 4, shift 157, rank 269, integer relations 0
```

The same queue record was then replayed with graph optimization enabled and a
180-second timeout:

```text
queue 0, Shape B, --cuso-graph on, 180s:
  timeout, multiplicity 9, shift 1016, no lattice relation stage reached
```

Graph optimization explores much deeper shift sets for this two-variable
partial-low600 model, but the run still did not reach a useful integer
relation, root, or factor within the bounded budget.

This is useful cost evidence, not a no-root certificate or factor.  The current
signal is that Shape B partial-low600 is wired correctly and the
no-intermediate option reaches larger lattices, while graph optimization reaches
larger shift sets.  Neither option toggle is enough at this fixed budget.

The first Shape C queue record was also replayed with graph-off and a
90-second timeout:

```text
queue 8, Shape C lower-origin:
  timeout, multiplicity 1, rank 768, integer relations 0

queue 8, Shape C --cuso-no-intermediate:
  timeout after building rank 4541 at multiplicity 1
```

Two fixed-budget D projections were collected with `PHASE_SEEDS=19` and one
iteration per shape:

```text
D_265_48_362_16: one 64-literal projection, variable_bits 506
D_265_64_362_8:  one 72-literal projection, variable_bits 498
```

Both D records were replayed with graph-off and a 60-second timeout:

```text
D_265_48_362_16: timeout, multiplicity 1, rank 1024, integer relations 0
D_265_64_362_8:  timeout, multiplicity 1, rank 1024, integer relations 0
```

Two adaptive fixed-budget E projections were collected with `PHASE_SEEDS=29`
and one iteration per shape:

```text
E_150_4_265_84_362_16: one 104-literal projection, variable_bits 466
E_150_4_265_84_362_32: one 120-literal projection, variable_bits 450
```

Both E records were replayed with graph-off and a 180-second timeout:

```text
E_150_4_265_84_362_16: timeout, multiplicity 3, shift 148, rank 192, integer relations 0
E_150_4_265_84_362_32: timeout, multiplicity 3, shift 243, rank 128, integer relations 0
```

Three narrower F projections were then collected to bracket the fully fixed
low600 leaf behavior:

```text
F_150_4_265_84_362_50: timeout at 120s, multiplicity 2, shift 259, rank 192, integer relations 0
F_150_4_265_84_362_54: timeout at 120s, multiplicity 1, rank 512, integer relations 0
F_150_4_265_84_362_58: soft_no_root in 1.991s, multiplicity 3, shift 5, rank 8, integer relations 8, roots 0
```

The F58 result is the expected leaf-oracle baseline: once all low600 holes are
fixed, cuso can cheaply derive integer relations and reject this random branch.
F50/F54 still time out, so leaving even 8 or 4 low bits unfixed keeps this
model outside the current bounded budget.

The same F58 leaf was then run inside the prefix-SAT loop with
`--learn-soft-no-root`.  A four-iteration smoke learned four 146-literal
blocking clauses from four cuso `soft_no_root` leaf checks.  After adding
JSONL logging for integrated cuso results, a two-iteration verification wrote:

```text
records: 2 prefix projections, 2 unique keys
result_status: soft_no_root=2
learned_clause events: 2
prefix build: 6.891s, 244,510 vars, 1,550,074 clauses
SAT solve per projection: 3.716s, 6.372s
cuso leaf checks: 2.352s, 1.552s
queue: logs/ct07_cuso/programmatic_f58_learn_loop_jsonl_20260706/projections.jsonl
results: logs/ct07_cuso/programmatic_f58_learn_loop_jsonl_20260706/results.jsonl
```

This proves the conservative solver-in-loop plumbing works for fully fixed
low600 leaves: SAT projection, cuso check, JSONL oracle result, learned clause,
and next projection.  It is not the desired broad-clause result because each
learned clause still fixes all 146 low600 unknown bits.

The next bridge is `leaf-certify`: for a parent shape with a small remaining
low variable, enumerate that small low variable and run the cheap F58 leaf
oracle for every completion.  This avoids asking cuso to solve the broader
multivariate parent shape directly.

For the earlier F54 projection that timed out as a two-variable cuso problem,
`leaf-certify` enumerated all 16 values of `u416_4`:

```text
shape: F_150_4_265_84_362_54
fixed bits: 142
leaf completions: 16
leaf statuses: soft_no_root=16
elapsed: 29.546s
result: leaf_exhausted_soft_no_root
results: logs/ct07_cuso/partial_f54_leaf_enum_20260706/results.jsonl
```

The same oracle was then wired into `ct07_programmatic_low600_sat_cas` as
`--oracle leaf-certify`.  A one-iteration F54 prefix-SAT loop produced a
142-literal projection, exhausted all 16 F58 leaves in 29.323s, and promoted
the parent projection to a learned clause:

```text
prefix build: 7.330s, 244,510 vars, 1,550,074 clauses
SAT solve: 4.447s
parent clause length: 142
leaf_status_counts: soft_no_root=16
learned_clause kind: leaf_exhausted_soft_no_root_promoted
queue: logs/ct07_cuso/programmatic_f54_leaf_loop_20260706/projections.jsonl
results: logs/ct07_cuso/programmatic_f54_leaf_loop_20260706/results.jsonl
```

This is the first practical path below the 146-literal leaf clause: F54 avoids
the 120-second direct cuso timeout by paying 16 cheap leaf checks.  It is still
classified as soft evidence unless the F58 leaf oracle is accepted as a
soundness gate, but it gives a concrete route for SAT-in-loop broad pruning.

The same idea was extended to the earlier F50 projection.  Direct F50 cuso had
timed out at 120s with two variables (`u412_8,z600`), but the new
`scripts/run_partial_leaf_shards.sh` runner split the 8-bit low variable into
eight parallel shards of 32 F58 leaf checks:

```text
shape: F_150_4_265_84_362_50
fixed bits: 138
leaf completions: 256
shards: 8 x 32
leaf_status: soft_no_root=256
summary_status: leaf_range_exhausted_soft_no_root=8
covered: 256
missing: 0
sum of per-leaf cuso elapsed: 1016.935s
per-shard elapsed range: 125.492s .. 128.816s
results: logs/ct07_cuso/partial_f50_leaf_shards_20260706/
```

This gives a second concrete widening step: a 138-literal F50 parent projection
can be exhausted by parallel leaf checks, even though the direct two-variable
cuso model did not finish inside the bounded run.  It is computationally much
heavier than F54, so this is best used as a selected-candidate verifier or a
parallel batch oracle, not as the default inner SAT callback yet.

The queue-level runner `scripts/run_programmatic_low600_queue_leaf_shards.sh`
now applies the same leaf-shard verifier directly to prefix-SAT projection
queues.  It computes `leaf_total = 2^(146-fixed_bits)` from each queue record
and skips records above `MAX_LEAF_BITS`.

It reproduced the first F50 result from the queue and also handled a second
fresh F50 projection collected with `phase_seed=89`:

```text
existing F50 queue record:
  queue_records=1
  leaf_records=256
  leaf_status=soft_no_root=256
  covered=256 missing=0
  logs=logs/ct07_cuso/programmatic_f50_leaf_queue_20260706/

new F50 queue record:
  key=F_150_4_265_84_362_50:63219137835160904486950493292694358870302
  queue_records=1
  leaf_records=256
  leaf_status=soft_no_root=256
  covered=256 missing=0
  logs=logs/ct07_cuso/programmatic_f50_leaf_queue_new_20260706/leaf_shards/
```

This makes selected F50 broad-pruning queue records repeatable without manually
copying fixed ids into a one-off command.

Verified leaf-shard parents can now be exported as reusable learned-clause
records with `experiments/export_leaf_learned.py`:

```bash
python3 experiments/export_leaf_learned.py \
  logs/ct07_cuso/programmatic_f50_leaf_queue_new_20260706/projections.jsonl \
  logs/ct07_cuso/programmatic_f50_leaf_queue_new_20260706/leaf_shards \
  --output-jsonl logs/ct07_cuso/learned_leaf_20260706/f50_new_learned.jsonl
```

The exported record has `record_type=learned_projection_clause`, the original
F50 fixed segments, and `source_status=all_leaves_soft_no_root`.  The prefix-SAT
runner accepts it with `--load-learned-jsonl`, and
`scripts/run_programmatic_low600_queue_collect.sh` forwards the same option via
`LOAD_LEARNED_JSONL`.

Using the same F50 `phase_seed=89` that originally produced fixed id
`63219137835160904486950493292694358870302`, loading
`f50_new_learned.jsonl` added one 138-literal clause and forced the next
projection to a different fixed id:

```text
loaded_learned_clauses: records=1, clauses=1, literals=138
old fixed id: 63219137835160904486950493292694358870302
new fixed id: 237443709698681397780198292297759683135774
queue: logs/ct07_cuso/learned_leaf_20260706/collector_load_smoke/projections.jsonl
```

This closes the loop for the current F50 workflow:

```text
prefix-SAT projection queue
-> sharded F58 leaf verification
-> learned_projection_clause JSONL
-> next prefix-SAT collection with learned clauses loaded
```

The loop can now be run as a single bounded cycle:

```bash
INITIAL_LEARNED_JSONL='logs/ct07_cuso/learned_leaf_20260706/f50_old_learned.jsonl logs/ct07_cuso/learned_leaf_20260706/f50_new_learned.jsonl' \
  SHAPES='F_150_4_265_84_362_50' PHASE_SEEDS='89' ITERATIONS=1 \
  QUEUE_LIMIT=1 MAX_LEAF_BITS=8 LEAF_SHARDS=8 LEAF_TIMEOUT=180s \
  bash scripts/run_programmatic_low600_learned_leaf_cycle.sh \
  logs/ct07_cuso/learned_leaf_cycle
```

The driver writes `learned_all.jsonl`, collects the next projection with all
previous learned clauses loaded, shard-verifies the selected queue records, and
merges any new `learned_projection_clause` records without duplicate keys.

One hundred thirty-two 2026-07-06 F50 parent projections have now been exported as learned
clauses.  Two came from the initial leaf-queue checks, four came from
continuing `phase_seed=89`, and one hundred twenty-six more used `phase_seed=97`,
`phase_seed=101`, `phase_seed=103`, `phase_seed=107`, `phase_seed=109`,
`phase_seed=113`, `phase_seed=127`, `phase_seed=131`, `phase_seed=149`, and
`phase_seed=151`, `phase_seed=157`, `phase_seed=163`, `phase_seed=167`, and
`phase_seed=173`, `phase_seed=179`, `phase_seed=181`, `phase_seed=191`, and
`phase_seed=193`, `phase_seed=197`, `phase_seed=199`, `phase_seed=211`, and
`phase_seed=223`, `phase_seed=227`, `phase_seed=229`, `phase_seed=233`, and
`phase_seed=239`, `phase_seed=241`, `phase_seed=251`, `phase_seed=257`, and
`phase_seed=263`, `phase_seed=269`, `phase_seed=271`, `phase_seed=277`, and
`phase_seed=281`, `phase_seed=283`, `phase_seed=293`, `phase_seed=307`, and
`phase_seed=311`, `phase_seed=313`, `phase_seed=317`, `phase_seed=331`, and
`phase_seed=337`, `phase_seed=347`, `phase_seed=349`, `phase_seed=353`, and
`phase_seed=359`, `phase_seed=367`, `phase_seed=373`, `phase_seed=379`, and
`phase_seed=383`, `phase_seed=389`, `phase_seed=397`, `phase_seed=401`, and
`phase_seed=409`, `phase_seed=419`, `phase_seed=421`, `phase_seed=431`, and
`phase_seed=433`, `phase_seed=439`, `phase_seed=443`, `phase_seed=449`, and
`phase_seed=457`, `phase_seed=461`, `phase_seed=463`, `phase_seed=467`, and
`phase_seed=479`, `phase_seed=487`, `phase_seed=491`, `phase_seed=499`, and
`phase_seed=503`, `phase_seed=509`, `phase_seed=521`, `phase_seed=523`, and
`phase_seed=541`, `phase_seed=547`, `phase_seed=557`, `phase_seed=563`, and
`phase_seed=569`, `phase_seed=571`, `phase_seed=577`, `phase_seed=587`, and
`phase_seed=593`, `phase_seed=599`, `phase_seed=601`, `phase_seed=607`, and
`phase_seed=613`, `phase_seed=617`, `phase_seed=619`, `phase_seed=631`, and
`phase_seed=641`, `phase_seed=643`, `phase_seed=647`, `phase_seed=653`, and
`phase_seed=659`, `phase_seed=661`, `phase_seed=673`, `phase_seed=677`, and
`phase_seed=683`, `phase_seed=691`, `phase_seed=701`, `phase_seed=709`, and
`phase_seed=719`, `phase_seed=727`, `phase_seed=733`, `phase_seed=739`, and
`phase_seed=743`, `phase_seed=751`, `phase_seed=757`, `phase_seed=761`, and
`phase_seed=769`, `phase_seed=773`, `phase_seed=787`, `phase_seed=797`, and
`phase_seed=809`, `phase_seed=811`, `phase_seed=821`, `phase_seed=823`, and
`phase_seed=827`, `phase_seed=829`, `phase_seed=839`, `phase_seed=853`, and
`phase_seed=857`, `phase_seed=859`, `phase_seed=863`, `phase_seed=877`, and
`phase_seed=881` to diversify the fixed
`150:4 + 265:84 + 362:50` projection:

```text
initial:
  F_150_4_265_84_362_50:206958246085563392169967157750982704324455
  F_150_4_265_84_362_50:63219137835160904486950493292694358870302
seed 89:
  F_150_4_265_84_362_50:237443709698681397780198292297759683135774
  F_150_4_265_84_362_50:150331423766921151133574392795227021003038
  F_150_4_265_84_362_50:324555995630441644426822191800292345268510
  F_150_4_265_84_362_50:19662994869280781163638543541428027803934
seed 97:
  F_150_4_265_84_362_50:105851071626254621238772737597800257631724
seed 101:
  F_150_4_265_84_362_50:172685572118155380379483846503171990793339
seed 103:
  F_150_4_265_84_362_50:47912146931635876251786862052849705335547
seed 107:
  F_150_4_265_84_362_50:171626801001643729947479492630250737748980
seed 109:
  F_150_4_265_84_362_50:185787797694225596835920538827514542912784
seed 113:
  F_150_4_265_84_362_50:268575149210774055687716376129961860550942
seed 127:
  F_150_4_265_84_362_50:47621549389091022516732349721778636948542
seed 131:
  F_150_4_265_84_362_50:161103982774720327522694061598817331111451
seed 149:
  F_150_4_265_84_362_50:136614063231845194878735791327733717469088
seed 151:
  F_150_4_265_84_362_50:159967903107378161740773502653254586356989
seed 157:
  F_150_4_265_84_362_50:243035995709807415508379355156831740134563
seed 163:
  F_150_4_265_84_362_50:101418305207644867403015404230536583480210
seed 167:
  F_150_4_265_84_362_50:101002021137977209848571018512278087366594
seed 173:
  F_150_4_265_84_362_50:336988953839654562779030293408100033075923
seed 179:
  F_150_4_265_84_362_50:332877831292966825909767063583000815561107
seed 181:
  F_150_4_265_84_362_50:207747695134385567137067548444225724092051
seed 191:
  F_150_4_265_84_362_50:136287977672327417572696559640355602250495
seed 193:
  F_150_4_265_84_362_50:69216211945789452138372421920171526848477
seed 197:
  F_150_4_265_84_362_50:42514618683059760587636618383515889124513
seed 199:
  F_150_4_265_84_362_50:225118146674086662553217941637980494859498
seed 211:
  F_150_4_265_84_362_50:291366194003502243338249968224307589648367
seed 223:
  F_150_4_265_84_362_50:235103430162809368600496259330818384459744
seed 227:
  F_150_4_265_84_362_50:310779354483176312406845414022912699286814
seed 229:
  F_150_4_265_84_362_50:328562852980269510720892543315087027526995
seed 233:
  F_150_4_265_84_362_50:331035948150805257207984693093201013576407
seed 239:
  F_150_4_265_84_362_50:9828903748261674167711523264401848951308
seed 241:
  F_150_4_265_84_362_50:119547639499996098587223164904868214929895
seed 251:
  F_150_4_265_84_362_50:161497458182288547776457461144110019379771
seed 257:
  F_150_4_265_84_362_50:39806792078506792582676781138501444679280
seed 263:
  F_150_4_265_84_362_50:58225196035802260108031472556653139221759
seed 269:
  F_150_4_265_84_362_50:261658997914226806573185258882192705252211
seed 271:
  F_150_4_265_84_362_50:279561568023937787558549851087318192592147
seed 277:
  F_150_4_265_84_362_50:347941732889304968995181230589780541019685
seed 281:
  F_150_4_265_84_362_50:294392260686605079505447602560229228603125
seed 283:
  F_150_4_265_84_362_50:244234297524121252901820106819600993951161
seed 293:
  F_150_4_265_84_362_50:283166135416663424685999558053711437429813
seed 307:
  F_150_4_265_84_362_50:269864209830697358899835422214130982119515
seed 311:
  F_150_4_265_84_362_50:163157992651198461294073657007633626452777
seed 313:
  F_150_4_265_84_362_50:257069972341111889103654963408205859922637
seed 317:
  F_150_4_265_84_362_50:20734479030173473588460030135374012520377
seed 331:
  F_150_4_265_84_362_50:70210082985451063322977557023459156813259
seed 337:
  F_150_4_265_84_362_50:55609359150129825328322879722875244612012
seed 347:
  F_150_4_265_84_362_50:191219357737083707442209505002298232931610
seed 349:
  F_150_4_265_84_362_50:138176150822200845171836996440195986495826
seed 353:
  F_150_4_265_84_362_50:256964765256280964651991501878717136607832
seed 359:
  F_150_4_265_84_362_50:138234705218244304990689250995989748447974
seed 367:
  F_150_4_265_84_362_50:66342330645386172232656380315439529043348
seed 373:
  F_150_4_265_84_362_50:18267778511849206573868870596883995634876
seed 379:
  F_150_4_265_84_362_50:99370142501295987638547265211691455862749
seed 383:
  F_150_4_265_84_362_50:289589336742687309518172213979620746569479
seed 389:
  F_150_4_265_84_362_50:330375196834013443439518004716905735264439
seed 397:
  F_150_4_265_84_362_50:107372668278824248913468787472259870107312
seed 401:
  F_150_4_265_84_362_50:242716839659512548682863512457156493775832
seed 409:
  F_150_4_265_84_362_50:94772997809772241986601648433422070905191
seed 419:
  F_150_4_265_84_362_50:111423533070305824179194900392032283912999
seed 421:
  F_150_4_265_84_362_50:10271275775724581862009432901163588073655
seed 431:
  F_150_4_265_84_362_50:6182965717744050020424985968138785180277
seed 433:
  F_150_4_265_84_362_50:98119457876866267597707901161416290385218
seed 439:
  F_150_4_265_84_362_50:109342809900005586306026583100068576079200
seed 443:
  F_150_4_265_84_362_50:87717259576459061637240131165676896239673
seed 449:
  F_150_4_265_84_362_50:171329509146416712962539753584846633345992
seed 457:
  F_150_4_265_84_362_50:321240077014516676553208767035628552292245
seed 461:
  F_150_4_265_84_362_50:308915179868457171952518507993943001650548
seed 463:
  F_150_4_265_84_362_50:82540719090718841040274191386404144500154
seed 467:
  F_150_4_265_84_362_50:103850001621923794277439968758682616023872
seed 479:
  F_150_4_265_84_362_50:166572119262559262196135822618046478153047
seed 487:
  F_150_4_265_84_362_50:256306309353868783049597448280816320337855
seed 491:
  F_150_4_265_84_362_50:227992833189040348523014498207825679247346
seed 499:
  F_150_4_265_84_362_50:202255872634798640903036926467957608553036
seed 503:
  F_150_4_265_84_362_50:126221484222362086995836695482668141925129
seed 509:
  F_150_4_265_84_362_50:31288209670933674871783625234479067690151
seed 521:
  F_150_4_265_84_362_50:54547521455678124213229097248166259808196
seed 523:
  F_150_4_265_84_362_50:194083862051337577541527266541025809922159
seed 541:
  F_150_4_265_84_362_50:279044928058204952751073334876044712585376
seed 547:
  F_150_4_265_84_362_50:187626469713829910168018139979169220338724
seed 557:
  F_150_4_265_84_362_50:165979226512333201034344964653862365351859
seed 563:
  F_150_4_265_84_362_50:346602493347316161489535031273048944416559
seed 569:
  F_150_4_265_84_362_50:46184162344841842370515463151690376107535
seed 571:
  F_150_4_265_84_362_50:249702658610698455454292123037178388396340
seed 577:
  F_150_4_265_84_362_50:144124490817238003174877745399113200667265
seed 587:
  F_150_4_265_84_362_50:259994902810392370298847518216144667264861
seed 593:
  F_150_4_265_84_362_50:276589177008743234027116497542668338707890
seed 599:
  F_150_4_265_84_362_50:82706523984672554718372227428016563431761
seed 601:
  F_150_4_265_84_362_50:223452146304939014331140966472496212398366
seed 607:
  F_150_4_265_84_362_50:193492907073990771131716185247182339759881
seed 613:
  F_150_4_265_84_362_50:205612981396645559612340250732658860471476
seed 617:
  F_150_4_265_84_362_50:118921675350876686268006971836781207429447
seed 619:
  F_150_4_265_84_362_50:113528743472603213503075699427424632084299
seed 631:
  F_150_4_265_84_362_50:197653177217086472204891930758204177250121
seed 641:
  F_150_4_265_84_362_50:269817496579995673176696392816306462333386
seed 643:
  F_150_4_265_84_362_50:152943060709586614062883888155133578740059
seed 647:
  F_150_4_265_84_362_50:286472472276780477865835083510498076173410
seed 653:
  F_150_4_265_84_362_50:267751414635428834954295292986593512324947
seed 659:
  F_150_4_265_84_362_50:320473315729406861502345729472311652009571
seed 661:
  F_150_4_265_84_362_50:158008089355856976969433426316462545799347
seed 673:
  F_150_4_265_84_362_50:8427709294261803386321788735236569220803
seed 677:
  F_150_4_265_84_362_50:236917633059510998527023233690921836632213
seed 683:
  F_150_4_265_84_362_50:100297902546302573984628620438338700528014
seed 691:
  F_150_4_265_84_362_50:4347800797050250136736166978613548520865
seed 701:
  F_150_4_265_84_362_50:49883885323696751871570309489601323521761
seed 709:
  F_150_4_265_84_362_50:42267924169500081837549610098567063913475
seed 719:
  F_150_4_265_84_362_50:43506964899764974287530006435447456061270
seed 727:
  F_150_4_265_84_362_50:76346617189170551604890383479396803885527
seed 733:
  F_150_4_265_84_362_50:190607154015708906027283003898727822419814
seed 739:
  F_150_4_265_84_362_50:233728209177449675173742611650684764651686
seed 743:
  F_150_4_265_84_362_50:161798904558952805123285887244925739443739
seed 751:
  F_150_4_265_84_362_50:125550498044850481810450910627033503477431
seed 757:
  F_150_4_265_84_362_50:183730408421053665133175055062690343202114
seed 761:
  F_150_4_265_84_362_50:181015242849926199912146872729872335934942
seed 769:
  F_150_4_265_84_362_50:7308113842936127491440655913437383831149
seed 773:
  F_150_4_265_84_362_50:192854747761032256269225327805099798713748
seed 787:
  F_150_4_265_84_362_50:217692254255113515519741030565894850952772
seed 797:
  F_150_4_265_84_362_50:95344023410978861892316270512069340420775
seed 809:
  F_150_4_265_84_362_50:109605123399930722399677316261755048107762
seed 811:
  F_150_4_265_84_362_50:37284145606735128329191698893225217682749
seed 821:
  F_150_4_265_84_362_50:239666158254263825370440478270225910326439
seed 823:
  F_150_4_265_84_362_50:139468241690907306083533268416161778143142
seed 827:
  F_150_4_265_84_362_50:138588977315886063113751657217628088703826
seed 829:
  F_150_4_265_84_362_50:126777694876206874288264420329496332066275
seed 839:
  F_150_4_265_84_362_50:311389213346072636102540149405412656239915
seed 853:
  F_150_4_265_84_362_50:337884940331690486575032892627232380267831
seed 857:
  F_150_4_265_84_362_50:259198180902820899380558416088827484846317
seed 859:
  F_150_4_265_84_362_50:290899975945046294611724265500365692818659
seed 863:
  F_150_4_265_84_362_50:194557369437289671468227425640387778468566
seed 877:
  F_150_4_265_84_362_50:155077476951113628069237846553896237989078
seed 881:
  F_150_4_265_84_362_50:90888320214480259957373632578859682516921
leaf_status: soft_no_root=256
summary_status: leaf_range_exhausted_soft_no_root=8
latest learned_all: 132 records
logs:
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed89_iter3/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed89_iter4/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed89_iter5_6/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed97_101/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed103_107/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed109_113/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed127_131/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed149_151/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed157_163/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed167_173/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed179_181/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed191_193/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed197_199/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed211_223/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed227_239/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed241_263/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed269_281/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed283_311/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed313_337/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed347_359/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed367_383/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed389_409/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed419_433/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed439_457/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed461_479/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed487_503/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed509_541/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed547_569/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed571_593/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed599_613/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed617_641/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed643_659/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed661_683/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed691_719/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed727_743/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed751_769/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed773_809/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed811_827/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed829_857/
  logs/ct07_cuso/learned_leaf_cycle_20260706_seed859_881/
latest seed859/881 cycle:
  queue_records=4
  leaf_records=1024
  summary_status: leaf_range_exhausted_soft_no_root=32
  export: learned=4, skipped=0
  merged total=132
```

So far the partial-low600 evidence favors single-variable leaf checks and
focus-template soft scoring over broader multivariate cuso clauses.  Shape B
remains the most useful broad-clause baseline to keep wired into
`ct07_programmatic_low600_sat_cas`; Shape C/D/E/F are coverage and cost
experiments, not sound no-root certificates or factors.

Set `CUSO_TIMEOUT=<seconds>` to bound each Sage+cuso child process:

```bash
CUSO_TIMEOUT=180 bash scripts/smoke_cuso_modes.sh logs/ct07_cuso/smoke_180s
CUSO_TIMEOUT=180 bash scripts/smoke_partial_low600_cuso.sh logs/ct07_cuso/partial_180s
CUSO_TIMEOUT=180 bash scripts/smoke_mixed_shape_cuso.sh logs/ct07_cuso/mixed_180s
CUSO_TIMEOUT=30s bash scripts/run_cuso_origin_matrix.sh logs/ct07_cuso/origin_matrix_30s
CUSO_TIMEOUT=60s bash scripts/run_cuso_option_sweep.sh logs/ct07_cuso/options_60s
CUSO_TIMEOUT=30s PLANTED_BITS=128 bash scripts/smoke_planted_cuso.sh logs/ct07_cuso/planted_128
PLANTED_BITS_LIST='64 96 128' CUSO_TIMEOUT=30s \
  bash scripts/run_planted_cuso_scale_sweep.sh logs/ct07_cuso/planted_scale_sweep
LOW_FREE_WIDTHS=0,1,2,4 CUSO_TIMEOUT=15s \
  bash scripts/run_planted_lowfree_sweep.sh logs/ct07_cuso/planted_lowfree_64
LOW_FREE_WIDTHS=1 CUSO_TIMEOUT=20s \
  bash scripts/run_planted_lowfree_option_sweep.sh logs/ct07_cuso/planted_lowfree1_options
LEAF_TOTAL=256 LEAF_SHARDS=8 LEAF_TIMEOUT=180s \
  bash scripts/run_partial_leaf_shards.sh \
  F_150_4_265_84_362_50 \
  206958246085563392169967157750982704324455 \
  logs/ct07_cuso/partial_f50_leaf_shards
QUEUE_LIMIT=1 MAX_LEAF_BITS=8 LEAF_SHARDS=8 LEAF_TIMEOUT=180s \
  bash scripts/run_programmatic_low600_queue_leaf_shards.sh \
  logs/ct07_cuso/programmatic_f50_smoke_20260706/projections.jsonl \
  logs/ct07_cuso/programmatic_f50_leaf_queue
python3 experiments/export_leaf_learned.py \
  logs/ct07_cuso/programmatic_f50_leaf_queue_new_20260706/projections.jsonl \
  logs/ct07_cuso/programmatic_f50_leaf_queue_new_20260706/leaf_shards \
  --output-jsonl logs/ct07_cuso/learned_leaf_20260706/f50_new_learned.jsonl
LOAD_LEARNED_JSONL=logs/ct07_cuso/learned_leaf_20260706/f50_new_learned.jsonl \
  SHAPES='F_150_4_265_84_362_50' PHASE_SEEDS='89' ITERATIONS=1 \
  bash scripts/run_programmatic_low600_queue_collect.sh \
  logs/ct07_cuso/learned_leaf_20260706/collector_load_smoke
python3 experiments/planted_focus_group_hm.py \
  --mode profile --prime-bits 64 --low-free-width 1 --m 4 --t 2 --rows 12
python3 experiments/planted_focus_group_hm.py \
  --mode sweep --prime-bits 96 --low-free-width 1 --m-values 10-12 --t-values 3-4 \
  --construction x-limited --x-shift-limit 1
python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep --prime-bits 96 --seeds 20260705-20260714 \
  --low-free-width 1 --m-values 11 --t-values 3 \
  --construction x-limited --x-shift-limit 1
python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep --prime-bits 96 --seeds 20260705-20260714 \
  --low-free-width 1 --m-values 11 --t-values 3 \
  --construction x-limited --x-shift-limit 1 --x-origin upper
python3 experiments/planted_focus_group_hm.py \
  --mode seed-sweep --prime-bits 96 --seeds 20260705-20260714 \
  --low-free-width 1 --m-values 11 --t-values 3 \
  --construction x-limited --x-shift-limit 1 \
  --template-families '0:0:3,0:1:3,1:0:2,1:1:0:2,1:1:2,2:0:1,2:1:0:1,2:1:1,2:1:1:1,2:2:1,3:0:0,3:1:0,3:1:0:0'
python3 experiments/planted_focus_group_hm.py \
  --mode sweep --prime-bits 64 --low-free-width 1 --m-values 12 --t-values 3-4 \
  --construction total-degree
python3 experiments/planted_focus_group_hm.py \
  --mode prune --prime-bits 64 --low-free-width 1 --m 12 --t 3
python3 experiments/planted_focus_group_hm.py \
  --mode prune --prime-bits 64 --low-free-width 1 --m 12 --t 3 \
  --min-support-coeff-bits 80
python3 experiments/planted_focus_group_hm.py \
  --mode drop-sweep --prime-bits 64 --low-free-width 1 --m 12 --t 3 \
  --drop-axes power,nscale
bash scripts/run_planted_focus_group_sweep.sh logs/ct07_cuso/planted_focus_group_64
PLANTED_BITS_LIST='64 80 96' \
  bash scripts/run_planted_focus_group_scale_sweep.sh logs/ct07_cuso/planted_focus_scale
CUSO_TIMEOUT=30s bash scripts/run_cuso_smoke_matrix.sh logs/ct07_cuso/matrix_30s
```

The provided shell runners export `DOT_SAGE=$PWD/.sage-tmp` and prepend
`$PWD/.sage-site` to `PYTHONPATH`, so they can see the local cuso install.
For manual `sage -python ...` commands, use the same environment:

```bash
DOT_SAGE="$PWD/.sage-tmp" PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}" \
  sage -python src/solve7_main.py --mode cuso --a 0 --b 1 --cuso-log INFO
```

On this machine, grouped cid0 with graph optimization entered cuso normally and
built a rank-275 dual lattice before a 180-second smoke timeout.  That confirms
the local Sage/cuso/flatter/msolve path is wired, but it is not a completed
factor search.

The 30-second matrix smoke gave these early signals:

```text
grouped graph off: rank 70, 0 integer relations
grouped graph on:  multiplicity 5, 261 shift polynomials
partial B graph on: multiplicity 6, 411 shift polynomials
partial B graph off: rank 75, 0 integer relations
partial C graph off: rank 512, 0 integer relations
partial A graph on: multiplicity 2, 3211 shift polynomials
```

The main `solve7_main.py` upper-low origin smoke also enters cuso cleanly:

```text
grouped upper-low graph off: rank 70, 0 integer relations
grouped upper-low graph on:  multiplicity 5, 261 shift polynomials
split upper-low graph on:    multiplicity 1
```

Treat Shape B as the first partial-low600 continuation candidate.  Shape A is
lower priority unless a different fixed-budget or cuso option controls the
shift blowup.

The 120-second graph-off follow-up reached:

```text
grouped cid0 graph off: rank 128, 0 integer relations
partial Shape B graph off: rank 128, 0 integer relations
```

A 300-second Shape B graph-off follow-up timed out after multiplicity 4,
157 shift polynomials, and rank 128 with 0 integer relations.  Shape B remains
the best broad-clause target conceptually, but this cost curve makes the
fixed-budget D shapes the next smoke targets before spending full sweeps on
grouped cid ranges.

A 180-second `D_265_64_362_8` graph-off smoke also timed out, with
multiplicity 1, rank 1536, and 0 integer relations.  Its lower variable mass did
not compensate for the jump to four cuso variables in this setting, so prefer
2-variable broad-clause shapes or mixed grouped shapes before expanding more
low600 holes into separate variables.

A 15-second mixed-shape graph-off smoke confirmed the new S0-S4 runner enters
cuso for each shape.  No shape produced integer relations in that short window:

```text
S0 grouped 2-var:              multiplicity 2, shift 21, rank 48
S1 exact 5-var:                multiplicity 1, rank 96
S2 low exact + high grouped:   multiplicity 1, rank 96
S3 low grouped + high exact:   multiplicity 1, rank 96
S4 low exact + high mixed:     multiplicity 1, rank 96
```

This keeps S0 and S2 as the most practical follow-up shapes for longer
bounded graph-off runs.  S3/S4 remain useful comparison points, but they start
with four variables and should not be expanded into full 256-edge sweeps first.

The first 120-second S2 graph-off follow-up timed out at multiplicity 1,
rank 256, and 0 integer relations.  That makes S2 useful as a shape-search
data point, but not yet a better default than S0/grouped for broad sweeps.

A 10-second option sweep over S0, S2, and partial Shape B found no integer
relations.  The strongest early signals were still only cost-shape signals:

```text
S0 graph off no-intermediate:        rank 70, 0 integer relations
S0 graph on:                         multiplicity 4, shift 151
S2 graph off no-intermediate:        rank 268, timed out before relation check
S2 graph on:                         multiplicity 2, shift 268
partial B graph off no-intermediate: rank 75, 0 integer relations
partial B graph on:                  multiplicity 5, shift 269
```

The S0 grouped no-intermediate path was rerun for 180 seconds on edge id 0.  It
confirmed 0 integer relations at ranks 21, 70, and 151, then generated the
multiplicity-4 ideal and built a rank-261 lattice before ending incomplete:

```text
S0 grouped 2-var, --cuso-no-intermediate, 180s:
  incomplete, multiplicity 4, shift 151, rank 261, integer relations 0
```

The S2 no-intermediate path was then rerun for 180 seconds on edge id 0 with
the repo-local Sage environment.  It completed the first rank-268 lattice
reduction and confirmed 0 integer relations, then generated the multiplicity-2
ideal and built a rank-1701 lattice before the run ended incomplete:

```text
S2 low-exact/high-grouped, --cuso-no-intermediate, 180s:
  incomplete, multiplicity 2, shift 268, rank 1701, integer relations 0
```

That is a useful cost boundary: the no-intermediate S2 model can pass the first
relation check, but the next lattice is too large for the current bounded run.

No option in that sweep made S0/S2/partial B look ready for a full range run.
Before spending time on full 256-edge sweeps, run
`scripts/run_cuso_option_sweep.sh` with a bounded timeout.  If no option creates
integer relations for S0/S2/partial B, the next useful work is focus-group or
planted-instance analysis rather than another blind cuso range sweep.

The planted true-branch check is in `experiments/planted_cuso_smoke.py`.  It
scales the challenge block geometry to a smaller prime size, constructs a known
RSA instance, and checks whether cuso recovers the planted factor for S0, S2,
partial-B, and a low-tail univariate baseline:

```bash
python3 experiments/planted_cuso_smoke.py --mode self-test --prime-bits 128
python3 experiments/planted_cuso_smoke.py --mode list-shapes --prime-bits 128

DOT_SAGE="$PWD/.sage-tmp" PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}" \
  sage -python experiments/planted_cuso_smoke.py \
  --mode cuso --shape S0_grouped_2var --prime-bits 128 \
  --cuso-log INFO --cuso-graph off --json
```

Use this as a true-branch retention smoke only.  A planted success does not make
real challenge no-root results sound; it only catches shape/API choices that
already fail on the known small case.

Initial 128-bit and 96-bit planted S0/S2/partial-B runs timed out without
recovering the planted factor, so the current planted lane first needs a smaller
64-bit low-tail univariate baseline before it can certify broader shapes.
That baseline succeeded on this machine:

```text
64-bit low_tail_univariate: factor, roots 1, elapsed 2.655s
64-bit S0 grouped 2-var:    timeout, multiplicity 3, shift 77, rank 128, intrel 0
96-bit S0/S2/partial-B:     timeout, no factor
128-bit S0/S2/partial-B:    timeout, no factor
```

A 64-bit S0 option sweep also failed to recover the planted factor:

```text
S0 graph on:                  timeout, multiplicity 6, shift 434
S0 graph off no-intermediate: timeout, multiplicity 3, shift 77, rank 164, intrel 0
S0 graph off allow-partial:   timeout, multiplicity 3, shift 77, no factor
```

A targeted 64-bit low-free-width sweep starts from the successful
`low_tail_univariate` model and opens only a suffix of the scaled low grouped
block as a second variable.  The first smoke on this machine was:

```text
low_free_width 0: factor, roots 1, elapsed 3.732s, rank 8, intrel 7
low_free_width 1: timeout, multiplicity 1, rank 256, intrel 0
low_free_width 2: timeout, multiplicity 1, rank 114, intrel 0
low_free_width 4: timeout, multiplicity 2, shift 39, rank 48, intrel 0
```

This is a sharper retention boundary than the earlier S0/S2 smoke: the
one-variable high-tail model is healthy, but adding even one low variable is
not retained within the current short cuso budget.  Treat partial-low600
no-root outputs as heuristic only until this low-free-width lane or a
focus-group lattice variant recovers planted true branches with low variables.

A 20-second option sweep for the first failing case, `low_free_width=1`, did
not find a rescuing option:

```text
graph off:                  timeout, multiplicity 1, rank 512, intrel 0
graph off allow-partial:    timeout, multiplicity 1, rank 768, intrel 0
graph on:                   timeout, multiplicity 1, no relation stage reached
graph auto:                 timeout, multiplicity 1, no relation stage reached
graph off no-intermediate:  timeout, multiplicity 1, no relation stage reached
```

So the current evidence points away from simply toggling cuso options.  The
next useful work is either a focus-group/fpylll planted lattice that can show
which two-variable shifts are useful, or an adaptive fixed-budget oracle that
keeps more low bits fixed before asking cuso for broad clauses.

The first focus-group HM profiler is in `experiments/planted_focus_group_hm.py`.
It builds a small HM-style two-variable lattice for the planted
`low_free_width=1` model, runs fpylll LLL with the transformation matrix, and
reports which input shift rows contribute to each reduced row plus whether
that row vanishes at the known true root.  The default construction is the
original y-only shift set `f^k*y^j`; `--construction x-limited` opens only
`0 <= i <= --x-shift-limit` in `f^k*x^i*y^j`, and `--construction
total-degree` opens all shifts with `i+j <= m-k`.  Initial checks:

```text
m=3,t=1: dim 10x10, LLL 0.002s, first 8 rows vanish=False
m=4,t=2: dim 15x15, LLL 0.012s, first 12 rows vanish=False
```

The `scripts/run_planted_focus_group_sweep.sh` runner sweeps the same planted
model over `m,t` and records the first reduced row that vanishes at the true
root.  On this machine, `m<=11,t<=6` produced no vanishing row.  The first
useful rows appeared at:

```text
m=12,t=3: dim 91x91, LLL 5.251s, vanish 1, first row 41, norm 198.5, support 85
m=12,t=4: dim 91x91, LLL 7.675s, vanish 1, first row 10, norm 263.0, support 81
```

The first support-family sketch for `m=12,t=3` includes high and low shift
families such as `0:0:3` through `0:9:3` plus `10:0:0`, `10:1:0`, `10:2:0`,
`11:0:0`, `11:1:0`, and `12:0:0`.  This gives the first concrete
focus-group pruning target: keep the shift families that contribute to these
vanishing rows and test a row-pruned planted lattice before trying the same
idea on the full challenge.

That support-family pruning preserves the planted relation:

```text
m=12,t=3: source 91x91, retained 85 families, pruned 85x91, vanish 1, row 41
m=12,t=4: source 91x91, retained 81 families, pruned 81x91, vanish 1, row 10
```

This is only a small row reduction, but it proves the focus-group measurement
is actionable: LLL can still recover the planted relation when non-supporting
input shifts are removed.  The next useful pruning step is to minimize those
support families further, ideally by dropping families or bands and checking
whether the same planted vanishing row remains.

A coefficient-threshold prune is much more fragile:

```text
m=12,t=3, min coeff bits 64: retained 85, vanish 1
m=12,t=3, min coeff bits 80: retained 84, vanish 0
m=12,t=4, min coeff bits 88: retained 81, vanish 1
m=12,t=4, min coeff bits 96: retained 80, vanish 0
```

So the first useful relation depends on nearly the full support set.  Further
focus-group reduction should drop structured bands or rebuild the shift set
around these families, not blindly discard small transformation coefficients.

Structured band drops confirm the same fragility.  For both `m=12,t=3` and
`m=12,t=4`, dropping any single `power` band or any single `n_scale_power` band
from the source support set lost all planted vanishing rows.  Even dropping the
single `power=12` family in the `m=12,t=3` case changed `85x91 -> 84x90` and
lost the relation.  That means the first useful planted relation is dense
across the HM shift families; the next focus-group step should not be naive
band deletion, but a new shift-set design that starts near the observed
support shape and changes the lattice construction itself.

A newer 96-bit seed sweep gives a better row-pruned target.  With
`--construction x-limited --x-shift-limit 1`, `m=11,t=3`, and seeds
`20260705..20260714`, the full `144x78` basis found planted vanishing rows in
5/10 lower-origin cases and 6/10 upper-origin cases.  The recurring support
set was the same 13-family template:

```text
0:0:3, 0:1:3,
1:0:2, 1:1:0:2, 1:1:2,
2:0:1, 2:1:0:1, 2:1:1, 2:1:1:1, 2:2:1,
3:0:0, 3:1:0, 3:1:0:0
```

Retaining only those 13 families reduces the lattice to `13x15`.  At 96 bits,
the template preserved 4/10 lower-origin planted successes and all 6/10
upper-origin planted successes.  The lower/upper success sets are disjoint and
cover all ten seeds.  The same template also scales to 128 bits: lower-origin
recovers 7/10 seeds, upper-origin recovers the remaining 3/10, and the union
again covers all seeds `20260705..20260714`.  LLL time drops to milliseconds
per seed.  This is the first strong focus-group template: it is far smaller
than the source basis, reusable across planted scales, and origin-sensitive in
the expected direction.  It is still only planted evidence, not a challenge
solver.

The same 13-family template is now wired to actual challenge E-shape
partial-low600 projections in `experiments/challenge_focus_template_hm.py`.
For the two SAT-derived E projections collected in
`logs/ct07_cuso/programmatic_queue_collect_e_20260706/projections.jsonl`, the
runner builds the expected `13x15` basis for both lower and upper origins:

```bash
python3 experiments/challenge_focus_template_hm.py --mode self-test
python3 experiments/challenge_focus_template_hm.py \
  --mode dry-run \
  --queue-jsonl logs/ct07_cuso/programmatic_queue_collect_e_20260706/projections.jsonl \
  --queue-limit 2 --x-origin both --json
```

The first real challenge smoke used `m=11,t=3`, the 13-family template, both
origins, and finite-field resultant/Groebner recovery with 56 primes:

```text
E_150_4_265_84_362_16 lower: LLL 0.113s, recover 41.335s, candidates 0
E_150_4_265_84_362_16 upper: LLL 0.152s, recover 30.605s, candidates 0
E_150_4_265_84_362_32 lower: LLL 0.081s, recover 30.538s, candidates 0
E_150_4_265_84_362_32 upper: LLL 0.069s, recover 34.231s, candidates 0
```

A short verbose recovery on the E32 lower-origin record found 0 modular root
pairs for the first primes in both resultant and Groebner paths.  This makes
the template runner useful as a fast heuristic exclusion/cost probe, but still
not a sound SAT learned-clause oracle.

The runner now stops early when all enabled recovery paths have 0 root pairs
for the same prime.  The previous full 56-prime scan is still available with
`--no-stop-on-empty-prime`, but the default fast path is much cheaper on false
projections.  On the same two E projections:

```text
E16 lower: empty prime 257, 2 path checks, recover 0.886s, candidates 0
E16 upper: empty prime 257, 2 path checks, recover 0.853s, candidates 0
E32 lower: empty prime 257, 2 path checks, recover 1.168s, candidates 0
E32 upper: empty prime 257, 2 path checks, recover 0.860s, candidates 0
```

A broader random-phase queue with seeds `31 37 41 43`, two iterations per seed,
and shapes E16/E32 collected 16 unique projections.  Running both origins for
each projection produced 32/32 fast heuristic exclusions at prime 257:

```text
records: 32 origin-runs
status: not_found 32
empty_prime: 257 for all 32
candidate_count total: 0
average recovery time: 0.983s
max recovery time: 2.329s
```

This gives a practical outer-loop filter for prioritizing SAT projections.

The same challenge focus-template runner was replayed on the first eight Shape
B projections from
`logs/ct07_cuso/programmatic_queue_collect_broad_20260706/projections.jsonl`.
Shape B has the two variables expected by the current runner:

```text
variables: 265:84, 600:424
basis: 13x15
records: 16 origin-runs
status: not_found 16
empty_prime: 257 for all 16
candidate_count total: 0
average recovery time: 0.618s
max recovery time: 0.898s
```

This is a cheap soft filter result, not a hard SAT clause.  The runner also
supports `--brute-small-vars WIDTH`, which enumerates small variables before
applying the same two-variable template.  This makes Shape C available by
brute-forcing `150:4` and keeping `362:58,600:424` as the two lattice
variables.  Replaying all eight Shape C projections from the broad queue,
queue indexes 8..15, produced:

```text
brute: 150:4 values 0..15
variables after brute: 362:58, 600:424
basis: 13x15
records: 256 origin-runs
queue coverage: 8 records x 16 brute values x 2 origins
status: not_found 256
empty_prime: 257 for all 256
candidate_count total: 0
average recovery time: 0.880s
median recovery time: 0.790s
max recovery time: 2.566s
```

It should still remain a soft score until a soundness gate exists for the
challenge-scale row-pruned HM polynomials.

The focus-template queue wrapper now forwards `FOCUS_BRUTE_SMALL_VARS`, so the
same batch runner can replay Shape C by brute-forcing the 4-bit `150:4`
variable before applying the two-variable template:

```bash
LOAD_LEARNED_JSONL=logs/ct07_cuso/learned_leaf_cycle_20260706_seed127_131/learned_all.jsonl \
  SHAPES='B C' PHASE_SEEDS='137 139' ITERATIONS=2 \
  bash scripts/run_programmatic_low600_queue_collect.sh \
  logs/ct07_cuso/programmatic_queue_collect_bc_seed137_139_20260706

FOCUS_TIMEOUT=90s FOCUS_BRUTE_SMALL_VARS=4 FOCUS_X_ORIGIN=both \
  FOCUS_RECOVER=both \
  bash scripts/run_programmatic_low600_queue_focus_template.sh \
  logs/ct07_cuso/programmatic_queue_collect_bc_seed137_139_20260706/projections.jsonl \
  logs/ct07_cuso/programmatic_focus_template_bc_seed137_139_20260706
```

This collected eight new random-phase projections while loading the 14 learned
F50 leaf clauses, then replayed them through the focus-template soft filter:

```text
queue: B=4, C=4, unique keys=8
learned input: 14 F50 projection clauses
Shape B origin-runs: 8
Shape C origin-runs: 128 = 4 records x 16 brute values x 2 origins
status: not_found 136/136
empty_prime: 257 for all 136
candidate_count total: 0
average recovery time: 0.697s
max recovery time: 1.149s
```

No factor/plaintext was found.  This strengthens the focus-template lane as a
fast projection ranker and confirms the Shape C brute replay path, but it is
still not a sound learned-clause oracle.

The next planted retention check widened the low-free variable in the same
two-variable planted model at 128 bits.  It compared the 13-family `13x15`
template against the full `x-limited --x-shift-limit 1`, `m=11,t=3`, `144x78`
source lattice on seeds `20260705..20260714`:

```text
low_free_width=2:
  template lower/upper union: 6/10 seeds
  full lower/upper union:     6/10 seeds

low_free_width=4:
  template lower/upper union: 3/10 seeds
  full lower/upper union:     3/10 seeds

low_free_width=8:
  template lower/upper union: 0/10 seeds
  full lower/upper union:     0/10 seeds
```

For these widths the sparse template exactly matches the full source lattice's
success set, so the failure is not caused by row pruning.  The current
`m=11,t=3` x-limited construction itself loses true-branch retention as the
low variable widens.  Since the challenge Shape C focus-template path leaves a
58-bit low variable after brute-forcing `150:4`, its fast empty-prime results
must remain soft scoring only; using them as hard SAT clauses needs a stronger
construction or a different soundness gate.

A direct width-8 parameter smoke on seed `20260705` then tried the full
x-limited source lattice with `m=11..13`, `t=3..5`, and both lower/upper
origins.  No setting produced a planted vanishing row:

```text
lower origin: 0/9 parameter points, max LLL 95.390s at m=13,t=5
upper origin: 0/9 parameter points, max LLL 98.833s at m=13,t=5
```

So the next planted work should not be a blind small `m/t` increase for the
same shift family.  It needs a different construction or a new way to split
the low variable before Shape C-style focus-template exclusions can become a
sound pruning oracle.

The planted runner now also accepts `--low-free-offset`, so the low-free window
can be moved inside the scaled low group instead of always using the upper
suffix.  On the same 128-bit seed `20260705`, `low_free_width=8`, `m=11,t=3`,
and offsets `0..11`, every lower/upper-origin run still had 0 planted
vanishing rows:

```text
offsets tested: 12
origin runs: 24
successes: 0
dimension per run: 144x78
max LLL time: 15.705s
```

So the width-8 failure is not just a bad window position.  For Shape C-style
hard pruning, simply choosing another 8-bit subwindow of the low block is not
enough under the current x-limited HM construction.

The scale check is in `scripts/run_planted_focus_group_scale_sweep.sh`.
Initial results:

```text
64-bit seed 20260705, low_free=1:
  m=12,t=3 and m=12,t=4 both produce one planted vanishing row.

80-bit seed 20260705/20260706, low_free=0:
  m=12,t=3/4 produce sparse rows, but this is a degenerate low-free branch.

80-bit seed 20260707, low_free=1:
  m=12,t=3 has no vanishing row.
  m=12,t=4 has one planted vanishing row and support-family pruning retains it.

96-bit seed 20260705, low_free=1:
  m=12,t=3/4/5/6 produce no planted vanishing rows.
  m=13,t=4/5 also produce no planted vanishing rows.
  m=14,t=4/5/6 and m=15,t=5/6 also produce no planted vanishing rows.
```

So the current `m=12` template scales from 64 to 80 bits only at the stronger
`t=4` setting, and does not reach 96 bits.  The first 96-bit `m=13` smoke also
misses, and `m=14/15` did not recover it either.  The next planted focus-group
work should change the shift-set construction before trying to project this
lattice shape onto the full 1024-bit challenge.

The first alternate shift-set check is `--construction total-degree`.  It is a
diagnostic option rather than a better default so far: after filtering zero
rows from the rectangular basis, 64-bit `m=12,t=3/4` finds one nonzero planted
vanishing row, but with `455x91` bases, 30.474s/105.379s LLL, and support
209/253.  The support-family prune for 64-bit `m=12,t=3` keeps the relation
after reducing `455x91 -> 209x91`.  On 80-bit seed `20260707`, total-degree
also needs `m=12,t=4` and takes 120.628s with support 253.  On 96-bit seed
`20260705`, total-degree `m=8..10,t=3..4` still finds no planted vanishing row.
So this confirms the construction is wired correctly, but not that it improves
the scale boundary over the y-only lattice.

The more useful alternate check is `--construction x-limited
--x-shift-limit 1`.  It adds a small number of `x` shifts without the full
total-degree row explosion.  On 64-bit seed `20260705`, it still first succeeds
only at `m=12,t=3/4`, with `169x91` bases, support 102/104, and support-family
pruning preserving the relation after `169x91 -> 102x90`.  On 80-bit seed
`20260707`, it again succeeds at `m=12,t=4`.  The important new signal is
96-bit: seed `20260705` succeeds at `m=11,t=3` with `144x78`, support 87, and
support-family pruning preserves the relation after `144x78 -> 87x77`.
Seed `20260707` is even sharper: `m=11,t=3` has seven planted vanishing rows,
and the first relation prunes from `144x78` to `13x15` while still vanishing.
Seed `20260706` still fails for `m=10..12,t=3..4`, so this is not a universal
template yet, but it is the first planted focus-group variant here that moves a
96-bit low-free instance past the y-only boundary.

The seed-sweep mode checks whether that 96-bit signal is stable:

```text
96-bit seeds 20260705..20260714, low_free=1, x-limited limit 1, m=11,t=3:
  success 5/10
  successes: 20260705, 20260707, 20260710, 20260711, 20260713
  failures:  20260706, 20260708, 20260709, 20260712, 20260714
  recurring families among successes include:
    0:0:3, 0:1:3, 1:0:2, 1:1:0:2, 1:1:2, 2:0:1, 2:1:0:1, 2:1:1
```

For those five failures, `x-limited --x-shift-limit 2` at `m=11,t=3`
recovered 0/5, and `x-limited --x-shift-limit 1` at `m=11,t=4` also recovered
0/5.  So the useful next search is not simply more `x` shifts or higher `t` at
the same `m`; it should compare nearby constructions around the recurring
small support families, especially the 13-family relation seen on seeds
`20260707`, `20260710`, `20260711`, and `20260713`.

`seed-sweep` also accepts `--template-families`, which builds only the named
shift families and then runs LLL.  The recurring 13-family sparse template

```text
0:0:3, 0:1:3,
1:0:2, 1:1:0:2, 1:1:2,
2:0:1, 2:1:0:1, 2:1:1, 2:1:1:1, 2:2:1,
3:0:0, 3:1:0, 3:1:0:0
```

uses a `13x15` lattice and reproduces 4/10 seeds:
`20260707`, `20260710`, `20260711`, and `20260713`.  It does not reproduce the
`20260705` dense relation.  A dense 87-family template extracted from the
`20260705` full x-limited relation uses an `87x77` lattice and reproduces all
five full x-limited successes while still rejecting the same five failures.
This means the current success/failure split is not caused by irrelevant extra
rows in the `144x78` source lattice.

The split is explained by the root's position in the one-bit `x` box.
`--x-origin upper` rewrites the first variable as `x' = X-1-x`, so a root at
the upper edge becomes a zero-near root.  On the same 96-bit seed range,
`x-limited --x-shift-limit 1 --x-origin upper`, `m=11,t=3` succeeds on the six
seeds where the original low-free bit is `1`:

```text
upper-origin successes:
  20260705, 20260706, 20260708, 20260709, 20260712, 20260714

zero-origin sparse-template successes:
  20260707, 20260710, 20260711, 20260713
```

The same 13-family sparse template at `13x15` reproduces the upper-origin
successes.  Together, zero-origin plus upper-origin sparse templates explain
all 10 planted 96-bit seeds.  This is the strongest planted focus-group result
so far: variable origin matters more than adding rows, and full challenge
fallbacks should try lower-edge and upper-edge origins for small/open low
variables before treating no-root as evidence.

This separates cuso/reconstruction correctness from broad-shape difficulty:
the known-low one-variable model works, while the grouped and partial broad
models still need focus-group or parameter analysis before their no-root
results can be trusted as pruning evidence.

After the environment smoke succeeds, run the partial-low600 oracle, then use
it inside `ct07_programmatic_low600_sat_cas` only after no-root soundness is
established.  Run `ct07_cuso_mixed_shape_search` in parallel; use
`ct07_focus_group_hm` and `ct07_cocert_clause_minimization` for local fallback
lattice pruning and ledger-to-clause repair rather than spending more time on
blind q-gap ledger batches.

The 2026-07-06 mixed-shape lower-origin smoke used graph-off and a 60-second
timeout per shape for edge id 0:

```text
S0_grouped_2var:             timeout, multiplicity 3, shift 70, rank 96,  integer relations 0
S1_exact_5var:               timeout, multiplicity 1,           rank 192, integer relations 0
S2_low_exact_high_grouped:   timeout, multiplicity 1,           rank 256, integer relations 0
S3_low_grouped_high_exact:   timeout, multiplicity 1,           rank 192, integer relations 0
S4_low_exact_high_mixed:     timeout, multiplicity 1,           rank 256, integer relations 0
```

The next lower-origin graph-off check used the neighboring edge id 1.  S0 for
this edge is covered by the grouped per-candidate runner (`cid 1`: timeout,
multiplicity 3, shift 70, rank 64, integer relations 0).  The exact and mixed
variants were then run with 60-second timeouts:

```text
S1_exact_5var:               timeout, multiplicity 1, rank 192, integer relations 0
S2_low_exact_high_grouped:   timeout, multiplicity 1, rank 256, integer relations 0
S3_low_grouped_high_exact:   timeout, multiplicity 1, rank 192, integer relations 0
S4_low_exact_high_mixed:     timeout, multiplicity 1, rank 192, integer relations 0
```

This edge-1 comparison again favors the two-variable grouped shape for cheap
coverage.  It does not show a relation-producing advantage for exact or mixed
3/4/5-variable shapes under short graph-off budgets.

S0 remains the lightest mixed/grouped shape in this short matrix.  A separate
S0 `--cuso-no-intermediate` smoke with a 180-second timeout advanced to
multiplicity 4, shift 151, and rank 261.  It confirmed 0 integer relations
through the preceding rank-151 reduction, then ended incomplete while reducing
the rank-261 lattice.  This mirrors the partial-low600 evidence: fewer
variables are more tractable, while option toggles alone have not yet produced
relations.

S0 origin variants were then compared with graph-off and a 90-second timeout
for the same edge id 0:

```text
S0 upper-low: timeout, multiplicity 3, shift 70, rank 128, integer relations 0
S0 upper-all: timeout, multiplicity 3, shift 70, rank 128, integer relations 0
```

For this edge candidate, flipping the low grouped variable or both grouped
variables did not improve on the lower-origin path.  The planted focus-group
origin result still matters for fallback lattices, but the current cuso S0
edge-0 evidence does not show an origin-only breakthrough.

To split the grouped path manually, candidate id is `cid = 16*high + low`:

```bash
bash scripts/run_cuso_grouped_candidates.sh 0 8 logs/ct07_cuso/grouped_candidates_0_8

bash scripts/run_cuso_range.sh 0 64
bash scripts/run_cuso_range.sh 64 128
bash scripts/run_cuso_range.sh 128 192
bash scripts/run_cuso_range.sh 192 256
```

Prefer `run_cuso_grouped_candidates.sh` for bounded exploration because it runs
one cid per Sage process and applies `CUSO_TIMEOUT` to each candidate.  The
older `run_cuso_range.sh` and 8-way runner apply timeout to a whole range, so a
single difficult candidate can prevent later candidates in that range from
being tried.

The first per-candidate grouped smoke covered cids 0..23 on 2026-07-06 with
`CUSO_TIMEOUT=45s` and `CUSO_GRAPH=off`:

```text
cid 0: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 1: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 2: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 3: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 4: timeout, multiplicity 3, shift 70, rank 48, intrel 0
cid 5: timeout, multiplicity 2, shift 21, rank 70, intrel 0
cid 6: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 7: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 8: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 9: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 10: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 11: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 12: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 13: timeout, multiplicity 3, shift 70, rank 48, intrel 0
cid 14: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 15: timeout, multiplicity 3, shift 70, rank 48, intrel 0
cid 16: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 17: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 18: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 19: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 20: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 21: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 22: timeout, multiplicity 3, shift 70, rank 64, intrel 0
cid 23: timeout, multiplicity 3, shift 70, rank 64, intrel 0
```

No factor or candidate root was found.  This is useful coverage evidence that
the external Sage+cuso/flatter path is installed and callable, but the
lower-origin grouped model still needs larger per-candidate budgets, origin
variants, or a mixed/partial shape before interpreting absence of roots.

The first upper-low grouped replay covered cids 0..7 with the same timeout and
graph settings, using `CUSO_UPPER_LOW_VARIABLES=1` so `x` is upper-origin and
`y` remains lower-origin:

```text
cid 0: timeout, multiplicity 3, shift 70, final rank 64, intrel 0
cid 1: timeout, multiplicity 3, shift 70, final rank 64, intrel 0
cid 2: timeout, multiplicity 3, shift 70, final rank 64, intrel 0
cid 3: timeout, multiplicity 3, shift 70, final rank 12, intrel 0
cid 4: timeout, multiplicity 3, shift 70, final rank 32, intrel 0
cid 5: timeout, multiplicity 3, shift 70, final rank 64, intrel 0
cid 6: timeout, multiplicity 3, shift 70, final rank 64, intrel 0
cid 7: timeout, multiplicity 3, shift 70, final rank 64, intrel 0
```

This did not produce a factor, root, or integer relation.  It also did not
improve on the lower-origin `cid 0..7` smoke; `cid 5` actually lost the
distinct lower-origin `multiplicity 2, shift 21, rank 70` profile and joined
the common multiplicity-3 profile.

The next upper-low grouped chunk covered cids 8..15 with the same 45-second
graph-off budget:

```text
cid 8:  timeout, multiplicity 3, shift 70, final rank 48, intrel 0
cid 9:  timeout, multiplicity 3, shift 70, final rank 48, intrel 0
cid 10: timeout, multiplicity 3, shift 70, final rank 48, intrel 0
cid 11: timeout, multiplicity 3, shift 70, final rank 48, intrel 0
cid 12: timeout, multiplicity 3, shift 70, final rank 48, intrel 0
cid 13: timeout, multiplicity 3, shift 70, final rank 48, intrel 0
cid 14: timeout, multiplicity 2, shift 21, final rank 70, intrel 0
cid 15: timeout, multiplicity 2, shift 21, final rank 70, intrel 0
```

The short upper-low run for cids 14 and 15 showed the same
`multiplicity 2, shift 21, rank 70` profile that had looked distinctive in the
lower-origin cid 5 run, so cid 14 was rerun with a 180-second timeout:

```text
cid 14 upper-low 180s:
  timeout, multiplicity 3, shift 70, final rank 151, integer relations 0
```

This also produced no root or factor.  The grouped path is still useful for
candidate coverage, but the current evidence argues against spending long
budgets on upper-low grouped candidates solely because their short-run profile
looks different; mixed/partial-low600 shapes are the better next use of cuso
time.

The selected lower-origin `cid 5` profile was then rerun with a 180-second
per-candidate timeout:

```text
cid 5 lower-origin 180s:
  timeout, multiplicity 4, shift 151, final rank 32, integer relations 0
```

The longer run confirmed that `cid 5` can progress past the 45-second
`multiplicity 2, shift 21, rank 70` point into the multiplicity-4 search, but it
still produced no integer relation, root, or factor.  This weakens the case for
spending long budgets on grouped candidates solely because their short-run
profiles differ; selected long grouped runs should now be compared against
partial-low600 Shape B or other lower-variable-count broad-clause shapes.

The split path also accepts `--a/--b`, but those are edge-id ranges rather than
grouped candidate-id ranges:

```bash
sage -python src/solve7_main.py --mode cuso-split \
  --cuso-split-brute-small-edges --a 0 --b 32
```

The mixed-shape cuso search is in `experiments/mixed_shape_cuso.py`.  It uses
the same `cid = 16*high + low` edge-id convention as grouped cuso and reports
both total variable mass and known-gap bits absorbed by grouped variables:

```bash
python3 experiments/mixed_shape_cuso.py --mode self-test
python3 experiments/mixed_shape_cuso.py --mode list-shapes
python3 experiments/mixed_shape_cuso.py --mode dry-run \
  --shape S2_low_exact_high_grouped --a 0 --b 1
python3 experiments/mixed_shape_cuso.py --mode dry-run \
  --shape S2_low_exact_high_grouped --a 0 --b 1 --upper-low-variables

sage -python experiments/mixed_shape_cuso.py \
  --mode cuso --shape S2_low_exact_high_grouped \
  --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
sage -python experiments/mixed_shape_cuso.py \
  --mode cuso --shape S2_low_exact_high_grouped \
  --a 0 --b 1 --upper-low-variables --cuso-log INFO --cuso-graph off --json
sage -python experiments/mixed_shape_cuso.py \
  --mode cuso --shape S3_low_grouped_high_exact \
  --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
sage -python experiments/mixed_shape_cuso.py \
  --mode cuso --shape S4_low_exact_high_mixed \
  --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
```

The mixed-shape runner supports the same origin variants as partial-low600:
`--upper-low-variables` flips variables ending before bit 600,
`--upper-variable START:WIDTH` flips a specific model variable, and
`--upper-all-variables` flips every variable.  Use lower and upper low-origin
smokes before interpreting a mixed-shape no-root as shape evidence.

The partial-low600 oracle is in `experiments/low600_partial_cuso.py`.  It can be
inspected without Sage:

```bash
python3 experiments/low600_partial_cuso.py --mode self-test
python3 experiments/low600_partial_cuso.py --mode list-shapes
python3 experiments/low600_partial_cuso.py --mode dry-run --shape B --a 0 --b 1
python3 experiments/low600_partial_cuso.py \
  --mode dry-run --shape B --a 0 --b 1 --upper-low-variables
```

Run one cuso smoke attempt per shape on the external Sage+cuso machine:

```bash
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape B --a 0 --b 1 --cuso-log INFO --cuso-graph on --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape B --a 0 --b 1 --upper-low-variables \
  --cuso-log INFO --cuso-graph on --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape C --a 0 --b 1 --cuso-log INFO --cuso-graph on --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape D_265_64_362_8 --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape D_265_48_362_16 --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape E_150_4_265_84_362_16 \
  --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape E_150_4_265_84_362_32 \
  --a 0 --b 1 --cuso-log INFO --cuso-graph off --json
sage -python experiments/low600_partial_cuso.py \
  --mode cuso --shape A --a 0 --b 1 --cuso-log INFO --cuso-graph on --json
```

`--upper-low-variables` applies the planted focus-group origin insight to the
actual partial-low600 oracle: low variables are represented as
`u' = 2^width-1-u`, while the high tail `z600` remains lower-origin.  Use this
alongside the default lower-origin run before interpreting a no-root result as
evidence.  `--upper-variable START:WIDTH` can flip a specific model variable,
and `--upper-all-variables` also flips `z600` for diagnostic sweeps.

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
