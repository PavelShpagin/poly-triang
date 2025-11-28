#!/usr/bin/env bash
set -euo pipefail

# Benchmark script for comparing triangulation methods including our new reflex algorithm

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/generated"
SCRIPTS_DIR="${ROOT}/scripts"
PUB_DIR="${ROOT}/publication"

mkdir -p "${RESULTS_DIR}" "${PUB_DIR}/figures"

echo "=== Polygon Triangulation Methods Benchmark ==="
echo ""

# Generate polygons with specific sizes for method comparison
echo "[1/6] Generating test polygons..."
python3 "${SCRIPTS_DIR}/generate_polygons.py" --output "${POLY_DIR}" --sizes 10 50 100 500 1000 2000 5000 10000

echo "[2/6] Building all executables..."
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1
cmake --build "${BUILD_DIR}" --target earcut_cli reflex_cli -j4 >/dev/null

echo "[3/6] Running benchmarks..."
BENCH_CSV="${RESULTS_DIR}/methods_benchmark.csv"
echo "algorithm,polygon,num_vertices,num_reflex,time_ms" > "${BENCH_CSV}"

echo ""
echo "Running benchmarks on all polygons..."
echo ""

run_earcut() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/earcut_${name}.tri"
  if log="$("${BIN_DIR}/earcut_cli" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "earcut,${name},${num_vertices},0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-14s %-20s %6s verts  %12s ms\n" "Earcut" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-14s %-20s FAILED\n" "Earcut" "${name}"
  fi
}

run_reflex() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  local output="${RESULTS_DIR}/reflex_${name}.tri"
  if log="$("${BIN_DIR}/reflex_cli" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    num_reflex="$(echo "${log}" | sed -n 's/.*reflex_count=\([0-9]*\).*/\1/p')"
    echo "reflex,${name},${num_vertices},${num_reflex},${time_ms}" >> "${BENCH_CSV}"
    printf "  %-14s %-20s %6s verts  r=%4s  %10s ms\n" "Reflex" "${name}" "${num_vertices}" "${num_reflex}" "${time_ms}"
  else
    printf "  %-14s %-20s FAILED\n" "Reflex" "${name}"
  fi
}

run_earclip_naive() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  # Skip for large polygons
  if [ "${num_vertices}" -gt 5000 ]; then
    printf "  %-14s %-20s %6s verts  SKIPPED (too slow)\n" "EarClipNaive" "${name}" "${num_vertices}"
    return
  fi
  
  local output="${RESULTS_DIR}/earclip_naive_${name}.tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_earclip_naive.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "earclip_naive,${name},${num_vertices},0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-14s %-20s %6s verts  %12s ms\n" "EarClipNaive" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-14s %-20s FAILED\n" "EarClipNaive" "${name}"
  fi
}

run_garey() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  # Skip for large polygons
  if [ "${num_vertices}" -gt 5000 ]; then
    printf "  %-14s %-20s %6s verts  SKIPPED (too slow)\n" "Garey" "${name}" "${num_vertices}"
    return
  fi
  
  local output="${RESULTS_DIR}/garey_${name}.tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_garey.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "garey,${name},${num_vertices},0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-14s %-20s %6s verts  %12s ms\n" "Garey" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-14s %-20s FAILED\n" "Garey" "${name}"
  fi
}

run_hertel() {
  local poly="$1"
  local name="$2"
  local num_vertices="$3"
  
  # Skip for large polygons
  if [ "${num_vertices}" -gt 5000 ]; then
    printf "  %-14s %-20s %6s verts  SKIPPED (too slow)\n" "Hertel" "${name}" "${num_vertices}"
    return
  fi
  
  local output="${RESULTS_DIR}/hertel_${name}.tri"
  if log="$(python3 "${SCRIPTS_DIR}/run_hertel.py" --input "${poly}" --output "${output}" 2>&1)"; then
    time_ms="$(echo "${log}" | sed -n 's/.*time_ms=\([0-9.]*\).*/\1/p')"
    echo "hertel,${name},${num_vertices},0,${time_ms}" >> "${BENCH_CSV}"
    printf "  %-14s %-20s %6s verts  %12s ms\n" "Hertel" "${name}" "${num_vertices}" "${time_ms}"
  else
    printf "  %-14s %-20s FAILED\n" "Hertel" "${name}"
  fi
}

for poly in "${POLY_DIR}"/*.poly; do
  name="$(basename "${poly}" .poly)"
  num_vertices="$(head -n1 "${poly}")"
  
  run_earcut "${poly}" "${name}" "${num_vertices}"
  run_reflex "${poly}" "${name}" "${num_vertices}"
  run_earclip_naive "${poly}" "${name}" "${num_vertices}"
  run_garey "${poly}" "${name}" "${num_vertices}"
  run_hertel "${poly}" "${name}" "${num_vertices}"
  echo ""
done

echo "[4/6] Results stored in ${BENCH_CSV}"
echo ""

echo "=== Summary ==="
python3 -c "
import pandas as pd
df = pd.read_csv('${BENCH_CSV}')
print('Average time (ms) by algorithm:')
print(df.groupby('algorithm')['time_ms'].mean().round(4).to_string())
print('\nResults by polygon size:')
pivot = df.groupby(['num_vertices', 'algorithm'])['time_ms'].mean().unstack()
# Reorder columns
cols = ['earcut', 'reflex', 'earclip_naive', 'garey', 'hertel']
cols = [c for c in cols if c in pivot.columns]
print(pivot[cols].round(4).to_string())
"
echo ""

echo "[5/6] Generating visualizations..."
python3 "${SCRIPTS_DIR}/visualize_methods.py"
echo ""

echo "[6/6] Done!"
echo "  Results:  ${RESULTS_DIR}/methods_benchmark.csv"
echo "  Figures:  ${PUB_DIR}/figures/"

