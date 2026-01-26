# Instructions for Completing the Paper

## Current Status

The implementation and paper are nearly complete. The following has been done:

1. **Algorithm Implementation** (`build/bin/reflex_cli`):
   - Chain-based sweep (output-sensitive in the number of extrema \(k\)): `methods/include/reflex_chain_triangulate.hpp`
   - Edge-based linked sweep fallback for high-reflex inputs: `methods/include/reflex_fast_linked.hpp`
   - Hybrid strategy in the CLI: chain-based for \(r < n/8\), linked fallback otherwise (`methods/src/reflex_cli.cpp`)

2. **Paper** (`paper/submission.tex`):
   - Complete theoretical + experimental writeup
   - Auto-includes benchmark tables from `paper/generated/`

## What Remains

### Step 1: Build and Test

```bash
cd /home/pavel/dev/poly-triang

# Build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make reflex_cli polypartition_mono_cli seidel_cli polypartition_hm_cli -j4
cd ..

# Quick test
python3 test_fast.py
```

### Step 2: Run Paper Benchmarks (authoritative for `submission.tex`)

```bash
python3 paper/benchmark.py --sizes 100 500 1000 2000 5000 10000 --polygons-per-config 5
```

This will:
- Generate polygons internally (convex, dent, star, random)
- Run all algorithms with multiple instances per configuration
- Write:
  - `paper/results.txt`
  - `paper/generated/benchmark_table.tex`
  - `paper/generated/benchmark_full.tex`

### Step 3: Verify Paper Includes Tables

`paper/submission.tex` includes the generated tables via:
- `\\input{generated/benchmark_table}`
- `\\input{generated/benchmark_full}`

### Step 4: Verify Claims

The paper claims:
- O(n + r log r) complexity for chain-based approach
- 10-20x speedup on convex polygons vs Garey
- Competitive or faster on low-r polygons
- Falls back gracefully to O(n log n) for high-r polygons

Verify these claims match the benchmark results. If not, adjust the discussion accordingly.

### Step 5: Final Compilation

```bash
cd paper
pdflatex submission.tex
pdflatex submission.tex
```

## Key Files

| File | Purpose |
|------|---------|
| `methods/include/reflex_chain_triangulate.hpp` | Chain-based sweep implementation |
| `methods/include/reflex_fast_linked.hpp` | Linked sweep fallback implementation |
| `methods/src/reflex_cli.cpp` | CLI wrapper for benchmarking |
| `paper/submission.tex` | Main paper LaTeX source |
| `paper/generated/benchmark_table.tex` | Generated benchmark table (random polygons) |
| `paper/generated/benchmark_full.tex` | Generated full benchmark table (all types) |
| `paper/benchmark.py` | Benchmark runner script (paper) |
| `test_fast.py` | Quick correctness and performance test |

## Algorithm Summary

### Complexity Analysis

| Phase | Chain-based (r < n/8) | Edge-based (r >= n/8) |
|-------|----------------------|----------------------|
| Classify vertices | O(n) | O(n) |
| Build chains | O(n) | -- |
| Sort events | O(r log r) | O(n log n) |
| Process events | O(r log r) | O(n log n) |
| Face extraction | O(n) | O(n) |
| Triangulation | O(n) | O(n) |
| **Total** | **O(n + r log r)** | **O(n log n)** |

### Critical Implementation Details

1. **Binary search with on-demand evaluation**: Chains maintain relative x-order, so we can binary search without recomputing all positions. Each query computes only O(log r) x-values.

2. **Lazy chain advancement**: Each vertex is visited at most once across all chain advancements, giving O(n) total amortized cost.

3. **Hybrid threshold**: r < n/8 uses chain-based (our algorithm), r >= n/8 uses edge-based (Garey-style). This balances asymptotic benefit against constant factor overhead.

## Expected Benchmark Results

Based on the algorithm design, expect:

| Polygon Type | r/n ratio | Expected Winner |
|--------------|-----------|-----------------|
| Convex | 0 | **Ours** (10-20x faster) |
| Spiral | ~0 | **Ours** (only one that works) |
| Low-reflex | < 1/8 | **Ours** (2-5x faster) |
| Random | ~1/2 | Competitive (within 2x) |
| Star | 1/2 | Competitive (within 2x) |

## Publication Target

**Computational Geometry: Theory and Applications (CGTA)**

The paper presents:
1. Novel O(n + r log r) algorithm - first practical sub-O(n log n) deterministic method
2. Complete correctness proof with slab-based visibility invariants
3. Rigorous complexity analysis
4. Comprehensive experimental evaluation
5. Practical hybrid implementation recommendation

## Contact

Authors:
- Pavel Shpagin (pavelandrewshpagin@knu.ua)
- Vasyl Tereschenko

Affiliation: Taras Shevchenko National University of Kyiv

