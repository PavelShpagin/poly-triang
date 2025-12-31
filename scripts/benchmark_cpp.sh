#!/usr/bin/env bash
# Comprehensive C++ benchmark for polygon triangulation
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${ROOT}/build/bin"
POLY_DIR="${ROOT}/polygons/generated"
RESULTS_DIR="${ROOT}/results"

mkdir -p "${RESULTS_DIR}"

echo "=== C++ Polygon Triangulation Benchmark ==="
echo ""

# Output CSV
CSV="${RESULTS_DIR}/cpp_benchmark.csv"
echo "algorithm,polygon_type,n,r,time_ms" > "${CSV}"

# Benchmark function
run_bench() {
    local algo="$1"
    local poly="$2"
    local ptype="$3"
    local n="$4"
    
    local output="${RESULTS_DIR}/${algo}_${ptype}_${n}.tri"
    
    if [ -x "${BIN_DIR}/${algo}_cli" ]; then
        for run in 1 2 3; do
            if log="$("${BIN_DIR}/${algo}_cli" --input "${poly}" --output "${output}" 2>&1)"; then
                time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
                r="$(echo "${log}" | sed -n 's/.*reflex_count=\([0-9]*\).*/\1/p')"
                [ -z "$r" ] && r="0"
                echo "${algo},${ptype},${n},${r},${time_ms}" >> "${CSV}"
            fi
        done
    fi
}

# Test sizes
SIZES="100 500 1000 2000 5000 10000 20000 50000"

# Test polygon types
for ptype in convex random spiral star; do
    echo "Testing ${ptype} polygons..."
    for n in ${SIZES}; do
        poly="${POLY_DIR}/${ptype}_${n}.poly"
        if [ -f "${poly}" ]; then
            echo "  n=${n}..."
            run_bench "reflex" "${poly}" "${ptype}" "${n}"
            run_bench "polytri" "${poly}" "${ptype}" "${n}"
        fi
    done
done

echo ""
echo "Results saved to ${CSV}"

# Generate summary
echo ""
echo "=== Summary (average of 3 runs) ==="
python3 -c "
import pandas as pd
df = pd.read_csv('${CSV}')
print()
print('Average time (ms) by algorithm and polygon type at n=10000:')
pivot = df[df['n'] == 10000].pivot_table(values='time_ms', index='polygon_type', columns='algorithm', aggfunc='mean')
if not pivot.empty:
    print(pivot.round(3).to_string())
"

