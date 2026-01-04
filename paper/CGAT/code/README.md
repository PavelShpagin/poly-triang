# CGAT Submission Code Artifact

This folder contains the **code artifact** accompanying the CGAT submission in `paper/CGAT/submission.tex`.

## Contents

- `ours/`: our implementation (output-sensitive in \(k\))
- `baselines/`: baseline implementations vendored for comparison (as used by the repository)

## Build (ours)

From `paper/CGAT/code/ours/`:

```bash
g++ -O3 -std=c++17 reflex_cli.cpp -o reflex_cli
```

Run:

```bash
./reflex_cli --input path/to/polygon.poly --output out.tri --algo chain
```

Notes:
- The CLI reports `k_count` (number of local maxima) and `reflex_count` (reflex vertices) for convenience.
- The implementation uses a **hybrid** strategy: it runs the chain-based method when \(k\) is small and falls back to a linked edge-based method otherwise (this is what the experimental harness uses).


