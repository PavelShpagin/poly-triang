# CGAT Submission Code Artifact

This folder contains the **code artifact** accompanying the CGAT submission in `paper/CGAT/submission.tex`.

## Contents

- `ours/`: our implementation (output-sensitive in \(k\))
- `baselines/`: baseline implementations vendored for comparison (as used by the repository)
  - Note: `baselines/` contains two vendored codebases (`polypartition` and `Seidel`), but the benchmark includes **three** baselines:
    - **Garey**: `polypartition_mono_cli.cpp` (monotone decomposition triangulation, built from `baselines/polypartition`)
    - **Seidel**: `seidel_cli.c` (built from `baselines/Seidel`)
    - **Hertel--Mehlhorn**: `cgal_hm_cli.cpp` (CGAL `approx_convex_partition_2`, not vendored under `baselines/`)

## Build (ours)

The top-level `paper/CGAT/build.sh` builds all binaries used by the benchmark.
To build only our CLI from `paper/CGAT/code/ours/`:

```bash
g++ -O3 -std=c++17 reflex_cli.cpp -o reflex_cli
```

Run:

```bash
./reflex_cli --input path/to/polygon.poly --output out.tri --algo chain_only
```

Notes:
- The CLI reports `k_count` (number of local maxima) and `reflex_count` (reflex vertices) for convenience.
- `--algo chain_only` runs the chain-based method (no fallback). This is what the benchmark harness uses.
- `--algo chain` is an alias of `chain_only` (kept for backwards compatibility).
- Use `--algo linked` to run only the edge-based linked method.


