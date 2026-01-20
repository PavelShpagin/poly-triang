#!/usr/bin/env bash
set -euo pipefail

CGAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[1/5] Build binaries..."
chmod +x "${CGAT_DIR}/build.sh" "${CGAT_DIR}/tools/"*.py
"${CGAT_DIR}/build.sh"

echo "[2/5] Sanity-check correctness (no crossings on sampled diagonals)..."
python3 "${CGAT_DIR}/tools/check_diagonal_validity.py" --n 500 --seed 0
python3 "${CGAT_DIR}/tools/check_diagonal_validity.py" --n 500 --seed 1

echo "[3/5] Run benchmarks (paper sizes)..."
TYPES="${CGAT_TYPES:-convex,dent,random}"
SIZES="${CGAT_SIZES:-500,1000,2000,5000,10000}"
RUNS="${CGAT_RUNS:-5}"
TIMEOUT="${CGAT_TIMEOUT:-10}"
PIN_CPU="${CGAT_PIN_CPU:-0}"

echo "  types:   ${TYPES}"
echo "  sizes:   ${SIZES}"
echo "  runs:    ${RUNS}"
echo "  timeout: ${TIMEOUT}s"
echo "  pin-cpu: ${PIN_CPU}"

python3 "${CGAT_DIR}/tools/benchmark_cgat.py" \
  --types "${TYPES}" \
  --sizes "${SIZES}" \
  --runs "${RUNS}" \
  --timeout "${TIMEOUT}" \
  --pin-cpu "${PIN_CPU}" \
  --bin-dir "${CGAT_DIR}/bin" \
  --out-csv "${CGAT_DIR}/generated/benchmark_results.csv" \
  --out-raw "${CGAT_DIR}/generated/benchmark_results_raw.csv"

echo "[4/5] Update CGAT LaTeX tables from raw benchmark CSV..."
python3 "${CGAT_DIR}/tools/generate_cgat_tables.py" \
  --input "${CGAT_DIR}/generated/benchmark_results_raw.csv" \
  --outdir "${CGAT_DIR}/generated"

echo "[5/5] Build submission PDF..."
cd "${CGAT_DIR}"
pdflatex -interaction=nonstopmode -halt-on-error submission.tex
pdflatex -interaction=nonstopmode -halt-on-error submission.tex

echo "Done:"
echo "  - ${CGAT_DIR}/generated/benchmark_table.tex"
echo "  - ${CGAT_DIR}/generated/benchmark_full.tex"
echo "  - ${CGAT_DIR}/submission.pdf"

