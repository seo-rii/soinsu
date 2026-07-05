# Experiments manifest

This folder contains intermediate scripts and logs produced during the working session.  They are included for auditability and future continuation, but `src/solve7_main.py` is the recommended entry point.

Notable files:

- `hm_*`, `solve_hm_actual.py`, `solve_try.py`: Herrmann-May/Coppersmith lattice experiments.
- `sat_*`, `run_sat_*`, `z3_*`: SAT and bit-multiplication encodings.
- `root_*`: resultant/Groebner/finite-field CRT root recovery experiments.
- `constants.py`, `solve_bits.py`: problem parsing and bit-structure analysis helpers.
- `basis_m3.txt`, `basis_m4.txt`: saved lattice basis artifacts from early 5-variable tests.
- `solve7_compact_previous.py`: the earlier compact script handed off before this package.

Architecture-specific compiled binaries were intentionally not required for the main workflow.  Rebuild any C++ helpers from source if needed.
