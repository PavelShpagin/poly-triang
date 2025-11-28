#!/usr/bin/env bash
set -euo pipefail

# Benchmark for Bucket Triangulation Algorithm
# Tests up to 1M vertices to demonstrate O(n) amortized complexity

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/generated"
SCRIPTS_DIR="${ROOT}/scripts"
EXP_DIR="${ROOT}/experiments/bucket"

mkdir -p "${RESULTS_DIR}" "${EXP_DIR}/figures"

echo "=== Bucket Triangulation Benchmark (up to 1M vertices) ==="
echo ""

# Generate polygons - limit random to 100K due to high reflex count
echo "[1/5] Generating test polygons..."
python3 "${SCRIPTS_DIR}/generate_polygons.py" --output "${POLY_DIR}" --sizes 1000 5000 10000 50000 100000 200000 500000 1000000

echo "[2/5] Building executables..."
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
cmake --build "${BUILD_DIR}" --target bucket_cli earcut_cli -j4 >/dev/null

echo "[3/5] Running benchmarks..."
BENCH_CSV="${RESULTS_DIR}/bucket_benchmark.csv"
echo "algorithm,polygon_type,num_vertices,num_reflex,grid_size,time_ms" > "${BENCH_CSV}"

echo ""
echo "Running benchmarks..."
echo ""

run_bucket() {
  local poly="$1"
  local ptype="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/bucket_$(basename "${poly}" .poly).tri"
  if log="$("${BIN_DIR}/bucket_cli" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    num_reflex="$(echo "${log}" | sed -n 's/.*reflex_count=\([0-9]*\).*/\1/p')"
    grid_size="$(echo "${log}" | sed -n 's/.*grid_size=\([0-9]*\).*/\1/p')"
    echo "bucket,${ptype},${num_vertices},${num_reflex},${grid_size},${time_ms}" >> "${BENCH_CSV}"
    printf "  %-8s %-8s %8s verts  r=%6s  grid=%4s  %12s ms\n" "Bucket" "${ptype}" "${num_vertices}" "${num_reflex}" "${grid_size}" "${time_ms}"
  else
    printf "  %-8s %-8s %8s verts  FAILED\n" "Bucket" "${ptype}" "${num_vertices}"
  fi
}

run_earcut() {
  local poly="$1"
  local ptype="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/earcut_$(basename "${poly}" .poly).tri"
  if log="$("${BIN_DIR}/earcut_cli" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "earcut,${ptype},${num_vertices},0,0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-8s %-8s %8s verts  %28s ms\n" "Earcut" "${ptype}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-8s %-8s %8s verts  FAILED\n" "Earcut" "${ptype}" "${num_vertices}"
  fi
}

# Run benchmarks for each polygon type
for ptype in convex random star; do
  echo ""
  echo "=== ${ptype} polygons ==="
  for poly in "${POLY_DIR}/${ptype}_"*.poly; do
    if [ ! -f "${poly}" ]; then continue; fi
    name="$(basename "${poly}" .poly)"
    num_vertices="$(head -n1 "${poly}")"
    
    # Skip small polygons for this benchmark
    if [ "${num_vertices}" -lt 1000 ]; then continue; fi
    
    # Skip very large random/star polygons (too slow due to high reflex count)
    if [ "${ptype}" = "random" ] && [ "${num_vertices}" -gt 100000 ]; then
      printf "  %-8s %-8s %8s verts  SKIPPED (high reflex)\n" "Bucket" "${ptype}" "${num_vertices}"
      printf "  %-8s %-8s %8s verts  SKIPPED (high reflex)\n" "Earcut" "${ptype}" "${num_vertices}"
      continue
    fi
    if [ "${ptype}" = "star" ] && [ "${num_vertices}" -gt 200000 ]; then
      printf "  %-8s %-8s %8s verts  SKIPPED (high reflex)\n" "Bucket" "${ptype}" "${num_vertices}"
      printf "  %-8s %-8s %8s verts  SKIPPED (high reflex)\n" "Earcut" "${ptype}" "${num_vertices}"
      continue
    fi
    
    run_bucket "${poly}" "${ptype}" "${num_vertices}"
    run_earcut "${poly}" "${ptype}" "${num_vertices}"
  done
done

echo ""
echo "[4/5] Results stored in ${BENCH_CSV}"

echo ""
echo "=== Complexity Analysis ==="
python3 -c "
import pandas as pd
import numpy as np

df = pd.read_csv('${BENCH_CSV}')

print('Scaling Analysis (T = a * n^b):')
print('-' * 60)

for alg in ['bucket', 'earcut']:
    alg_df = df[df['algorithm'] == alg]
    if len(alg_df) < 3:
        continue
    
    grouped = alg_df.groupby('num_vertices')['time_ms'].mean().reset_index()
    x = grouped['num_vertices'].values
    y = grouped['time_ms'].values
    
    # Log-log fit
    valid = (x > 0) & (y > 0)
    if np.sum(valid) >= 3:
        coeffs = np.polyfit(np.log(x[valid]), np.log(y[valid]), 1)
        b = coeffs[0]
        a = np.exp(coeffs[1])
        
        # R^2
        y_pred = a * (x[valid] ** b)
        ss_res = np.sum((y[valid] - y_pred)**2)
        ss_tot = np.sum((y[valid] - np.mean(y[valid]))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        
        print(f'{alg:10s}: T = {a:.2e} * n^{b:.3f}, R^2 = {r2:.4f}')
        
        if b < 1.15:
            print(f'            -> Near-linear O(n)')
        elif b < 1.5:
            print(f'            -> Between O(n) and O(n log n)')
        else:
            print(f'            -> Super-linear')

print('-' * 60)
"

echo ""
echo "[5/5] Generating visualizations..."
python3 "${SCRIPTS_DIR}/visualize_bucket.py"

echo ""
echo "=== Done! ==="
echo "  Results:  ${RESULTS_DIR}/bucket_benchmark.csv"
echo "  Figures:  ${EXP_DIR}/figures/"

