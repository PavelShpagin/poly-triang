# CGAT Submission Code Artifact

This folder contains the **code artifact** accompanying the CGAT submission in `paper/CGAT/submission.tex`.

## Contents

- `ours/`: our implementation (output-sensitive in \(k\))
- `baselines/`: baseline implementations vendored for comparison (as used by the repository)

## Build (ours)

The top-level `paper/CGAT/build.sh` builds all binaries used by the benchmark.
To build only our CLI from `paper/CGAT/code/ours/`:

```bash
g++ -O3 -std=c++17 reflex_cli.cpp -o reflex_cli
```

Run:

```bash
./reflex_cli --input path/to/polygon.poly --output out.tri --algo chain
```

Notes:
- The CLI reports `k_count` (number of local maxima) and `reflex_count` (reflex vertices) for convenience.
- `--algo chain` uses a **hybrid** strategy: it runs the chain-based method when \(k < \max(8, n/8)\) and falls back to the edge-based linked method otherwise (this is what the experimental harness uses).
- Use `--algo chain_only` to run the chain-based method with **no fallback**.
- Use `--algo linked` to run only the edge-based linked method.


