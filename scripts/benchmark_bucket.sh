#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/generated"
SCRIPTS_DIR="${ROOT}/scripts"

mkdir -p "${RESULTS_DIR}"

echo "=== Bucket Triangulation Benchmark (up to 1M vertices) ==="
echo ""

echo "[1/5] Generating test polygons..."
python3 "${SCRIPTS_DIR}/generate_polygons.py" --output "${POLY_DIR}" --sizes 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000

echo "[2/5] Building executables..."
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "${BUILD_DIR}" --target bucket_cli reflex_cli earcut_cli >/dev/null

echo "[3/5] Running benchmarks..."
BENCH_CSV="${RESULTS_DIR}/bucket_benchmark.csv"
echo "algorithm,polygon,num_vertices,reflex_count,grid_size,time_ms" > "${BENCH_CSV}"

echo "Running benchmarks..."
echo ""

# Define algorithms and their CLIs
ALGORITHMS=("Bucket" "Reflex" "Earcut")
CLI_BINS=("bucket_cli" "reflex_cli" "earcut_cli")

# Define polygon types
POLY_TYPES=("convex" "random" "star")

for ptype in "${POLY_TYPES[@]}"; do
  echo ""
  echo "=== ${ptype} polygons ==="
  for poly in "${POLY_DIR}/${ptype}_"*.poly; do
    name="$(basename "${poly}" .poly)"
    num_vertices="$(head -n1 "${poly}")"
    
    # Skip very large random/star polygons ONLY for Earcut due to extreme slowness
    if [[ "${ptype}" == "random" || "${ptype}" == "star" ]]; then
        if [[ "${num_vertices}" -ge 500000 ]]; then
             printf "  %-10s %-10s %10s verts  SKIPPED (high reflex)\n" "Earcut" "${ptype}" "${num_vertices}"
             # We continue to run Reflex and Bucket!
        fi
    fi

    for i in "${!ALGORITHMS[@]}"; do
      alg="${ALGORITHMS[i]}"
      bin="${CLI_BINS[i]}"
      output="${RESULTS_DIR}/${alg}_${name}.tri"
      
      # Skip logic
      if [[ "${alg}" == "Earcut" && "${num_vertices}" -ge 500000 && ("${ptype}" == "random" || "${ptype}" == "star") ]]; then
          continue
      fi
      
      log_output=""
      if [[ "${alg}" == "Bucket" || "${alg}" == "Reflex" ]]; then
        log_output="$("${BIN_DIR}/${bin}" --input "${poly}" --output "${output}")"
        time_ms="$(echo "${log_output}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
        reflex_count="$(echo "${log_output}" | sed -n 's/.*reflex_count=\([0-9]*\).*/\1/p')"
        grid_size="$(echo "${log_output}" | sed -n 's/.*grid_size=\([0-9]*\).*/\1/p')"
        printf "  %-10s %-10s %10s verts  r=%5s  grid=%4s  %10s ms\n" "${alg}" "${ptype}" "${num_vertices}" "${reflex_count}" "${grid_size}" "${time_ms}"
        echo "${alg},${name},${num_vertices},${reflex_count},${grid_size},${time_ms}" >> "${BENCH_CSV}"
      else # Earcut
        log_output="$("${BIN_DIR}/${bin}" --input "${poly}" --output "${output}")"
        time_ms="$(echo "${log_output}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
        printf "  %-10s %-10s %10s verts                      %10s ms\n" "${alg}" "${ptype}" "${num_vertices}" "${time_ms}"
        echo "${alg},${name},${num_vertices},0,0,${time_ms}" >> "${BENCH_CSV}" # 0 for reflex_count and grid_size
      fi
    done
  done
done

echo ""
echo "[4/5] Results stored in ${BENCH_CSV}"
echo ""

echo "=== Complexity Analysis ==="
python3 -c "
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

df = pd.read_csv('${BENCH_CSV}')
print('Scaling Analysis (T = a * n^b):')
print('------------------------------------------------------------')

def power_law(x, a, b):
    return a * (x ** b)

def log_power_law(x, log_a, b):
    return log_a + b * np.log(x)

for alg in df['algorithm'].unique():
    alg_data = df[df['algorithm'] == alg].groupby('num_vertices')['time_ms'].mean().reset_index()
    alg_data = alg_data[alg_data['time_ms'] > 0]
    
    if len(alg_data) < 2:
        continue

    x_data = alg_data['num_vertices'].values
    y_data = alg_data['time_ms'].values
    
    try:
        popt, pcov = curve_fit(log_power_law, x_data, np.log(y_data), p0=[np.log(y_data[0]), 1])
        log_a_fit, b_fit = popt
        a_fit = np.exp(log_a_fit)
        
        y_pred = power_law(x_data, a_fit, b_fit)
        r_squared = r2_score(np.log(y_data), np.log(y_pred))
        
        print(f'{alg:<10}: T = {a_fit:.2e} * n^{b_fit:.3f}, R^2 = {r_squared:.4f}')
        if alg == 'bucket':
            print('            -> Super-linear')
        elif alg == 'earcut':
            print('            -> Between O(n) and O(n log n)')
        elif alg == 'reflex':
            print('            -> O(n + r log r) or O(n) for convex')
    except Exception as e:
        print(f'Error fitting power law for {alg}: {e}')
print('------------------------------------------------------------')
"

echo "[5/5] Generating visualizations..."
python3 "${SCRIPTS_DIR}/visualize_bucket.py"
echo ""

echo "=== Done! ==="
echo "  Results:  ${RESULTS_DIR}/bucket_benchmark.csv"
echo "  Figures:  ${ROOT}/experiments/bucket/figures/"
