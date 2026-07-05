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

# Install cuso into the Sage Python environment.
git clone https://github.com/keeganryan/cuso.git
cd cuso
sage -pip install .
```

Then run from this asset package root:

```bash
sage -python src/solve7_main.py --mode analyze
sage -python src/solve7_main.py --mode cuso --a 0 --b 256 | tee logs/cuso_full.log
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
