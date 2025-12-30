#!/usr/bin/env bash
set -euo pipefail

cd /home/pavel/dev/poly-triang

echo "=== COMPREHENSIVE BENCHMARK: Earcut vs Reflex ==="
echo ""

echo "--- CONVEX POLYGONS (r=0, expect O(n) for Reflex) ---"
for sz in 1000 5000 10000 50000 100000; do
    echo "n=$sz:"
    ./build/bin/earcut_cli -i /tmp/bench_poly/convex_$sz.poly -o /tmp/out.tri
    ./build/bin/reflex_cli -i /tmp/bench_poly/convex_$sz.poly -o /tmp/out.tri
done
echo ""

echo "--- STAR POLYGONS (high r, stress test) ---"
for sz in 1000 5000 10000 50000 100000; do
    echo "n=$sz:"
    ./build/bin/earcut_cli -i /tmp/bench_poly/star_$sz.poly -o /tmp/out.tri
    ./build/bin/reflex_cli -i /tmp/bench_poly/star_$sz.poly -o /tmp/out.tri
done
echo ""

echo "--- RANDOM POLYGONS (moderate r) ---"
for sz in 1000 5000 10000 50000 100000; do
    echo "n=$sz:"
    ./build/bin/earcut_cli -i /tmp/bench_poly/random_$sz.poly -o /tmp/out.tri
    ./build/bin/reflex_cli -i /tmp/bench_poly/random_$sz.poly -o /tmp/out.tri
done
echo ""

echo "=== DONE ==="

