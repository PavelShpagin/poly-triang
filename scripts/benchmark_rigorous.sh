#!/usr/bin/env bash
set -euo pipefail

# Rigorous benchmark for polygon triangulation methods
# Tests all methods on uniform sampling of polygon sizes

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/generated"
SCRIPTS_DIR="${ROOT}/scripts"
PUB_DIR="${ROOT}/publication"

mkdir -p "${RESULTS_DIR}" "${PUB_DIR}/figures"

echo "=== Rigorous Polygon Triangulation Benchmark ==="
echo ""

# Generate polygons with extended sizes including 100K
echo "[1/5] Generating test polygons (including 100K)..."
python3 "${SCRIPTS_DIR}/generate_polygons.py" --output "${POLY_DIR}" --sizes 100 500 1000 2000 5000 10000 20000 50000 100000

echo "[2/5] Building executables..."
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
cmake --build "${BUILD_DIR}" --target earcut_cli reflex_cli -j4 >/dev/null

echo "[3/5] Running benchmarks..."
BENCH_CSV="${RESULTS_DIR}/rigorous_benchmark.csv"
echo "algorithm,polygon_type,num_vertices,num_reflex,time_ms" > "${BENCH_CSV}"

echo ""
echo "Running benchmarks..."
echo ""

# Function to run earcut
run_earcut() {
  local poly="$1"
  local ptype="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/earcut_$(basename "${poly}" .poly).tri"
  if log="$("${BIN_DIR}/earcut_cli" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "earcut,${ptype},${num_vertices},0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-10s %-8s %7s verts  %12s ms\n" "Earcut" "${ptype}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-10s %-8s %7s verts  FAILED\n" "Earcut" "${ptype}" "${num_vertices}"
  fi
}

# Function to run reflex
run_reflex() {
  local poly="$1"
  local ptype="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/reflex_$(basename "${poly}" .poly).tri"
  if log="$("${BIN_DIR}/reflex_cli" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    num_reflex="$(echo "${log}" | sed -n 's/.*reflex_count=\([0-9]*\).*/\1/p')"
    echo "reflex,${ptype},${num_vertices},${num_reflex},${time_ms}" >> "${BENCH_CSV}"
    printf "  %-10s %-8s %7s verts  r=%6s  %10s ms\n" "Reflex" "${ptype}" "${num_vertices}" "${num_reflex}" "${time_ms}"
  else
    printf "  %-10s %-8s %7s verts  FAILED\n" "Reflex" "${ptype}" "${num_vertices}"
  fi
}

# Function to run naive ear clipping (Python) - skip for large polygons
run_naive() {
  local poly="$1"
  local ptype="$2"
  local num_vertices="$3"
  
  # Skip for polygons > 2000 vertices (too slow)
  if [ "${num_vertices}" -gt 2000 ]; then
    return
  fi
  
  local output="${RESULTS_DIR}/naive_$(basename "${poly}" .poly).tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_earclip_naive.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "naive,${ptype},${num_vertices},0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-10s %-8s %7s verts  %12s ms\n" "Naive" "${ptype}" "${num_vertices}" "${time_ms}"
  fi
}

# Run benchmarks for each polygon type
for ptype in convex random spiral star; do
  echo ""
  echo "=== ${ptype} polygons ==="
  for poly in "${POLY_DIR}/${ptype}_"*.poly; do
    if [ ! -f "${poly}" ]; then continue; fi
    name="$(basename "${poly}" .poly)"
    num_vertices="$(head -n1 "${poly}")"
    
    run_earcut "${poly}" "${ptype}" "${num_vertices}"
    run_reflex "${poly}" "${ptype}" "${num_vertices}"
    run_naive "${poly}" "${ptype}" "${num_vertices}"
  done
done

echo ""
echo "[4/5] Results stored in ${BENCH_CSV}"

echo ""
echo "=== Summary by Algorithm ==="
python3 -c "
import pandas as pd
import numpy as np
df = pd.read_csv('${BENCH_CSV}')
print(df.groupby('algorithm')['time_ms'].agg(['mean', 'std', 'min', 'max']).round(3).to_string())
"

echo ""
echo "=== Summary by Size (Earcut vs Reflex) ==="
python3 -c "
import pandas as pd
df = pd.read_csv('${BENCH_CSV}')
pivot = df[df['algorithm'].isin(['earcut', 'reflex'])].pivot_table(
    index='num_vertices', 
    columns='algorithm', 
    values='time_ms', 
    aggfunc='mean'
)
if 'earcut' in pivot.columns and 'reflex' in pivot.columns:
    pivot['speedup'] = pivot['earcut'] / pivot['reflex']
print(pivot.round(3).to_string())
"

echo ""
echo "[5/5] Generating visualizations..."
python3 "${SCRIPTS_DIR}/visualize_rigorous.py"

echo ""
echo "=== Done! ==="
echo "  Results:  ${RESULTS_DIR}/rigorous_benchmark.csv"
echo "  Figures:  ${PUB_DIR}/figures/"

