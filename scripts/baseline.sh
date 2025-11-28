#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/generated"
SCRIPTS_DIR="${ROOT}/scripts"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${ROOT}/paper/figures"

echo "=== Polygon Triangulation Benchmark Suite ==="
echo ""

echo "[1/5] Generating polygons..."
python3 "${SCRIPTS_DIR}/generate_polygons.py" --output "${POLY_DIR}"

echo "[2/5] Configuring and building C++ baselines..."
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
cmake --build "${BUILD_DIR}" --target earcut_cli >/dev/null 2>&1 || true

echo "[3/5] Running benchmarks..."
BENCH_CSV="${RESULTS_DIR}/benchmark_results.csv"
echo "algorithm,polygon,num_vertices,time_ms" > "${BENCH_CSV}"

run_earcut() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/earcut_${name}.tri"
  if [ -x "${BIN_DIR}/earcut_cli" ]; then
    if log="$("${BIN_DIR}/earcut_cli" --input "${poly}" --output "${output}" 2>&1)"; then
      time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
      echo "earcut,${name},${num_vertices},${time_ms}" >> "${BENCH_CSV}"
      printf "  %-12s %-20s %6s verts  %10s ms\n" "Earcut" "${name}" "${num_vertices}" "${time_ms}"
    fi
  fi
}

run_garey() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  # Skip for very large polygons (too slow in Python)
  if [ "${num_vertices}" -gt 10000 ]; then
    printf "  %-12s %-20s %6s verts  SKIPPED (too slow)\n" "Garey" "${name}" "${num_vertices}"
    return
  fi
  
  local output="${RESULTS_DIR}/garey_${name}.tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_garey.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "garey,${name},${num_vertices},${time_ms}" >> "${BENCH_CSV}"
    printf "  %-12s %-20s %6s verts  %10s ms\n" "Garey" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-12s %-20s FAILED\n" "Garey" "${name}"
  fi
}

run_hertel() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  # Skip for very large polygons (too slow)
  if [ "${num_vertices}" -gt 10000 ]; then
    printf "  %-12s %-20s %6s verts  SKIPPED (too slow)\n" "Hertel" "${name}" "${num_vertices}"
    return
  fi
  
  local output="${RESULTS_DIR}/hertel_${name}.tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_hertel.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "hertel,${name},${num_vertices},${time_ms}" >> "${BENCH_CSV}"
    printf "  %-12s %-20s %6s verts  %10s ms\n" "Hertel" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-12s %-20s FAILED\n" "Hertel" "${name}"
  fi
}

run_earclip_naive() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  # Skip for very large polygons (O(n^2) is too slow)
  if [ "${num_vertices}" -gt 5000 ]; then
    printf "  %-12s %-20s %6s verts  SKIPPED (O(n^2) too slow)\n" "EarClipNaive" "${name}" "${num_vertices}"
    return
  fi
  
  local output="${RESULTS_DIR}/earclip_naive_${name}.tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_earclip_naive.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "earclip_naive,${name},${num_vertices},${time_ms}" >> "${BENCH_CSV}"
    printf "  %-12s %-20s %6s verts  %10s ms\n" "EarClipNaive" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-12s %-20s FAILED\n" "EarClipNaive" "${name}"
  fi
}

echo ""
echo "Running benchmarks on all polygons..."
echo ""

for poly in "${POLY_DIR}"/*.poly; do
  name="$(basename "${poly}" .poly)"
  num_vertices="$(head -n1 "${poly}")"
  
  run_earcut "${poly}" "${name}" "${num_vertices}"
  run_earclip_naive "${poly}" "${name}" "${num_vertices}"
  run_garey "${poly}" "${name}" "${num_vertices}"
  run_hertel "${poly}" "${name}" "${num_vertices}"
  echo ""
done

echo "[4/5] Results stored in ${BENCH_CSV}"

# Show summary
echo ""
echo "=== Summary ==="
python3 -c "
import pandas as pd
df = pd.read_csv('${BENCH_CSV}')
print('Average time (ms) by algorithm:')
print(df.groupby('algorithm')['time_ms'].mean().sort_values().to_string())
print()
print('Results by polygon size:')
pivot = df.pivot_table(values='time_ms', index='num_vertices', columns='algorithm', aggfunc='mean')
print(pivot.round(4).to_string())
"

echo ""
echo "[5/5] Generating visualizations..."
python3 "${SCRIPTS_DIR}/visualize.py" || echo "Visualization generation had issues"

echo ""
echo "=== Done! ==="
echo "  Results: ${BENCH_CSV}"
echo "  Figures: ${ROOT}/paper/figures/"
echo "  Paper:   ${ROOT}/paper/paper.tex"
