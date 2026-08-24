# External setup notes

## Preferred environment

Use a Linux machine with enough RAM for lattice reduction and polynomial solving.  The recommended path is SageMath plus cuso/flatter/msolve.

The cuso EUROCRYPT 2025 artifact lists these dependencies:

- SageMath 9.8 or newer recommended
- flatter
- msolve
- cuso from `https://github.com/keeganryan/cuso.git`

A typical install flow is:

```bash
# System packages vary by distribution.  Install SageMath first.
sudo apt-get update
sudo apt-get install -y git build-essential python3-pip sagemath

# Install flatter and msolve according to your OS/package manager.
# If packages are unavailable, build them from their upstream repositories.

# Install cuso into a project-local Sage Python target.  This avoids system
# site-packages and works without sudo.
cd /path/to/rsa_partial_leak_assets
mkdir -p .sage-tmp .sage-site
DOT_SAGE="$PWD/.sage-tmp" sage -pip install --target "$PWD/.sage-site" \
  git+https://github.com/keeganryan/cuso.git
```

Then run from this asset package root:

```bash
DOT_SAGE="$PWD/.sage-tmp" PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}" \
  sage -python src/solve7_main.py --mode analyze
DOT_SAGE="$PWD/.sage-tmp" PYTHONPATH="$PWD/.sage-site${PYTHONPATH:+:$PYTHONPATH}" \
  sage -python src/solve7_main.py --mode cuso --a 0 --b 256 | tee logs/cuso_full.log
```

The bundled shell runners set these two environment variables automatically.
Set `CUSO_TIMEOUT=<seconds>` when smoke-running shapes so a single difficult
candidate cannot capture the machine indefinitely:

```bash
CUSO_TIMEOUT=180 bash scripts/smoke_cuso_modes.sh logs/smoke_180s
CUSO_TIMEOUT=180 bash scripts/smoke_mixed_shape_cuso.sh logs/mixed_180s
CUSO_TIMEOUT=60s bash scripts/run_cuso_option_sweep.sh logs/options_60s
CUSO_TIMEOUT=30s PLANTED_BITS_LIST='64 96 128' \
  bash scripts/run_planted_cuso_scale_sweep.sh logs/planted_scale_sweep
```

For the programmatic SAT+CAS scaffold, install PySAT into a project-local
Python target:

```bash
cd /path/to/rsa_partial_leak_assets
mkdir -p .py-site
python3 -m pip install --target .py-site 'python-sat[pblib,aiger]>=1.8.dev13'
PYTHONPATH="$PWD/.py-site${PYTHONPATH:+:$PYTHONPATH}" \
  python3 experiments/programmatic_low600_sat_cas.py --mode self-test
```

## Parallel candidate splitting

Candidate id is `cid = 16*high + low`.  For 8 workers:

```bash
for i in $(seq 0 7); do
  a=$((i*32)); b=$(((i+1)*32));
  bash scripts/run_cuso_range.sh "$a" "$b" &
done
wait
```

## Fallback-only Dockerfile

`Dockerfile.local` is provided only for the pure-Python/fpylll fallback.  It does not install Sage/cuso.

```bash
docker build -f env/Dockerfile.local -t solve7-local .
docker run --rm -it solve7-local python3 src/solve7_main.py --mode analyze
docker run --rm -it solve7-local python3 src/solve7_main.py --mode local --a 0 --b 1
```
