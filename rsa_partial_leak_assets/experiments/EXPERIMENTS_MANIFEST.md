# Experiments manifest

This folder contains intermediate scripts and logs produced during the working session.  They are included for auditability and future continuation, but `src/solve7_main.py` is the recommended entry point.

Notable files:

- `hm_*`, `solve_hm_actual.py`, `solve_try.py`: Herrmann-May/Coppersmith lattice experiments.
- `low600_partial_cuso.py`: partial-low600 cuso broad-clause oracle for
  `ct07_partial_low600_cuso_broad_clause`, including dry-run shape summaries,
  self-tests, Sage+cuso execution for Shape A/B/C and fixed-budget variants,
  adaptive fixed-budget E variants, narrow F/leaf-oracle variants, and
  lower/upper-origin variable variants for low600 variables.  Its
  `leaf-certify` mode enumerates small remaining low variables and runs the
  one-variable F58 leaf oracle for every completion.  It supports
  `--leaf-start/--leaf-stop` range sharding for wider parents such as F50.
- `mixed_shape_cuso.py`: grouped/exact/hybrid cuso shape runner for
  `ct07_cuso_mixed_shape_search`, using the grouped cuso edge-id convention and
  reporting absorbed known-gap bits for each shape.  It also supports
  lower/upper-origin variable variants for low-side grouped or exact variables.
- `scripts/run_cuso_origin_matrix.sh`: configurable lower/upper-origin smoke
  matrix for selected mixed-shape and partial-low600 cuso shapes, summarized by
  `summarize_cuso_logs.py`.
- `scripts/smoke_cuso_modes.sh`, `scripts/run_cuso_grouped_8way.sh`,
  `scripts/run_cuso_grouped_candidates.sh`, and
  `scripts/run_cuso_split_edges_8way.sh`: main `src/solve7_main.py` cuso
  runners.  `run_cuso_grouped_candidates.sh` runs one grouped candidate per
  Sage process and per-candidate timeout, which is useful when a range-level
  timeout would get stuck on the first candidate.  They accept
  `CUSO_UPPER_LOW_VARIABLES=1`,
  `CUSO_UPPER_ALL_VARIABLES=1`, and space-separated `CUSO_UPPER_VARIABLES`
  entries to exercise the same lower/upper-origin idea on grouped and split
  challenge models.
- `programmatic_low600_sat_cas.py`: PySAT outer-loop scaffold for
  `ct07_programmatic_low600_sat_cas`.  It reuses `sat_mul_cnf.build()`, extracts
  Shape B/C low600 projections from SAT models, and calls the partial-low600
  oracle.  Its `prefix-loop` mode builds `p*q == N mod 2^k`, which can sample
  Shape B/C projections without waiting for a full factor model.  It can append
  projection records to `--queue-jsonl` and replay them with
  `--mode replay-queue`.  Replay supports dry-run, cuso, and focus-template
  HM filter oracles; integrated prefix-loop cuso and leaf-certify results are
  written to `--results-jsonl`; `soft_no_root` and leaf-exhaustion learning
  remain behind the explicit `--learn-soft-no-root` flag.
- `export_leaf_learned.py`: converts fully covered leaf-shard queue results
  into `learned_projection_clause` JSONL records.  The programmatic SAT runner
  can load these with `--load-learned-jsonl`, and the queue collector forwards
  them via `LOAD_LEARNED_JSONL`.
- `merge_learned_clauses.py`: merges learned-clause JSONL files while
  de-duplicating projection keys.  This keeps iterative leaf-shard cycles from
  reloading the same parent projection more than once.
- `scripts/smoke_programmatic_low600_sat_cas.sh`: self-test plus prefix-SAT
  Shape B/C dry-run smoke.  Set `RUN_CUSO=1` to add a bounded Sage+cuso
  integration smoke for Shape B.  It accepts `PROJECTION_PHASE` and
  `PHASE_SEED` to exercise zero/one/alternate/random projection phase hints,
  writes `projections.jsonl`, and dry-run replays the collected queue.
- `scripts/run_programmatic_low600_queue_cuso.sh`: per-record cuso replay
  runner for a projection queue.  It runs each queue item under its own timeout,
  writes separate logs, appends `results.jsonl`, and summarizes the logs with
  queue index/key metadata.
- `scripts/run_programmatic_low600_queue_focus_template.sh`: per-record
  focus-template HM replay runner for a projection queue.  It uses the same
  `replay-queue` interface without Sage, writes one log per queue item, appends
  `results.jsonl`, and emits both queue and focus summaries for fast soft
  scoring.  Set `FOCUS_BRUTE_SMALL_VARS=4` for Shape C-style queues so the
  4-bit `150:4` variable is enumerated before the two-variable template runs.
- `scripts/run_partial_leaf_shards.sh`: parallel runner for `leaf-certify`
  ranges.  It is useful for selected F50-style parent projections where the
  remaining low variable is too wide for one sequential inner-loop call but can
  be exhausted by sharded F58 leaf checks.
- `scripts/run_programmatic_low600_queue_leaf_shards.sh`: queue-level wrapper
  around `run_partial_leaf_shards.sh`.  It reads prefix projection records,
  derives the leaf completion count from `fixed_bits`, and runs selected F50/F54
  records through the sharded leaf verifier.
- `scripts/run_programmatic_low600_queue_collect.sh`: batch collector for
  prefix-SAT projection queues across shapes, projection phase policies, phase
  seeds, and iteration counts.  Set `LOAD_LEARNED_JSONL` to load exported
  learned projection clauses before each SAT solve.
- `scripts/run_programmatic_low600_learned_leaf_cycle.sh`: one-cycle or
  multi-cycle driver for the current F50 broad-clause workflow.  It merges
  initial learned files, collects new prefix projections with those clauses
  loaded, runs selected queue records through leaf shards, exports newly
  verified learned clauses, and writes a de-duplicated `learned_all.jsonl`.
- `summarize_projection_queue.py`: reports projection queue size, unique keys,
  duplicate count, shape/phase/prefix/clause-length distribution, and optional
  oracle-result status counts.
- `planted_cuso_smoke.py`: downscaled planted RSA/cuso smoke for checking
  true-branch retention before treating heuristic no-root observations as
  pruning evidence.
- `planted_lowfree_sweep.py`: downscaled true-branch retention boundary sweep
  that starts from the high-tail-only model and opens a controlled suffix or
  `--low-free-offset` window of the low grouped block as an additional cuso
  variable.  The companion
  `scripts/run_planted_lowfree_option_sweep.sh` runner compares graph and
  intermediate-size options for the first failing low-free-width case.
- `planted_focus_group_hm.py`: fpylll/Herrmann-May focus-group profiler for
  planted low-free models.  It records LLL transformation-matrix support for
  each reduced row and checks whether those rows vanish at the known true root.
  It also has support-family, coefficient-threshold, and structured band-drop
  row-prune modes, plus `--construction x-limited` and `--construction
  total-degree` for `f^k*x^i*y^j` shift-set diagnostics.  `--x-origin upper`
  flips the first variable to `X-1-x`, which lets the same sparse low-edge
  template handle upper-edge roots.  `--low-free-offset` is passed through to
  the planted model builder for position-dependent low-window checks.  The companion
  `scripts/run_planted_focus_group_sweep.sh`
  runner sweeps `m,t`, writes the first vanishing-row candidates, and records
  pruned-lattice retention checks.  `scripts/run_planted_focus_group_scale_sweep.sh`
  repeats the same checks across planted prime sizes; both runners accept the
  `CONSTRUCTION` and `X_ORIGIN` environment variables.
  `scripts/run_planted_focus_group_seed_sweep.sh`
  scans multiple planted seeds for one construction/parameter point and reports
  success counts plus recurring support families.  The seed sweep also accepts
  `--template-families` or the runner's `TEMPLATE_FAMILIES` environment variable
  to test fixed row-pruned support templates directly.
- `challenge_focus_template_hm.py`: applies the planted 13-family
  row-pruned HM template to actual challenge partial-low600 two-variable
  projections.  It can also enumerate small variables first with
  `--brute-small-vars`, which lets Shape C brute-force `150:4` and then reuse
  the same two-variable lattice on `362:58,600:424`.  It supports lower/upper
  low-variable origins, queue replay, LLL-only cost checks, finite-field CRT
  recovery through the existing resultant/Groebner helpers, and an empty-prime
  fast path for obvious projection exclusions.  Outputs are heuristic evidence
  only, not sound no-root certificates.
- `summarize_cuso_logs.py`: extracts exit status, multiplicity, shift counts,
  graph subset sizes, lattice rank, and root counts from cuso smoke logs.
- `sat_*`, `run_sat_*`, `z3_*`: SAT and bit-multiplication encodings.
- `root_*`: resultant/Groebner/finite-field CRT root recovery experiments.
- `constants.py`, `solve_bits.py`: problem parsing and bit-structure analysis helpers.
- `basis_m3.txt`, `basis_m4.txt`: saved lattice basis artifacts from early 5-variable tests.
- `solve7_compact_previous.py`: the earlier compact script handed off before this package.

Architecture-specific compiled binaries were intentionally not required for the main workflow.  Rebuild any C++ helpers from source if needed.
