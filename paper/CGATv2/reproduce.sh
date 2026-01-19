#!/usr/bin/env bash
set -euo pipefail

CGAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[1/4] Build binaries..."
chmod +x "${CGAT_DIR}/build.sh" "${CGAT_DIR}/tools/"*.py
"${CGAT_DIR}/build.sh"

echo "[2/5] Sanity-check correctness (no crossings on sampled diagonals)..."
python3 "${CGAT_DIR}/check_diagonal_validity.py" --n 500 --seed 0
python3 "${CGAT_DIR}/check_diagonal_validity.py" --n 500 --seed 1

echo "[3/5] Run benchmarks (paper sizes)..."
python3 "${CGAT_DIR}/tools/benchmark_cgat.py" \
  --types convex,dent,random,star \
  --sizes 500,1000,2000,5000,10000 \
  --runs 5 \
  --timeout 10 \
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

