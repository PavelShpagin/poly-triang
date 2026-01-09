#!/bin/bash
# Build and run the benchmark for O(n + r log r) polygon triangulation paper
#
# Compares our algorithm against:
#   - Garey et al. (polypartition): O(n log n) monotone decomposition
#   - Hertel-Mehlhorn (polypartition): O(n log n)
#   - Seidel: O(n log* n) randomized trapezoidal decomposition
#
# Usage:
#   ./run_benchmark.sh                 # default: sizes 1000,10000,100000, 10 runs
#   ./run_benchmark.sh --sizes 100000,500000 --runs 5

set -e

cd "$(dirname "$0")"

echo "============================================================"
echo "O(n + r log r) Polygon Triangulation Benchmark"
echo "============================================================"
echo ""

echo "[1/3] Building project..."
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. >/dev/null
make -j$(nproc) reflex_cli polypartition_mono_cli polypartition_hm_cli seidel_cli 2>&1 | grep -E "(Building|Linking|Built)" || true
cd ..
echo "      Build complete!"
echo ""

echo "[2/3] Running benchmark..."
python3 scripts/benchmark_paper.py "$@"

echo ""
echo "[3/3] Results saved to results/benchmark_results.csv"
echo ""
echo "============================================================"
echo "To update paper tables, manually transfer results to:"
echo "  paper/reflex_triang.tex (tab:results, tab:extended)"
echo "============================================================"
