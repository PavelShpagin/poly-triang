#!/usr/bin/env bash
set -euo pipefail

CGAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[1/4] Build binaries..."
chmod +x "${CGAT_DIR}/build.sh" "${CGAT_DIR}/tools/"*.py
"${CGAT_DIR}/build.sh"

echo "[2/4] Run benchmarks (paper sizes)..."
python3 "${CGAT_DIR}/tools/benchmark_cgat.py" \
  --types convex,dent,random,star \
  --sizes 100,500,1000,2000,5000,10000 \
  --runs 5 \
  --timeout 10 \
  --bin-dir "${CGAT_DIR}/bin" \
  --out-csv "${CGAT_DIR}/generated/benchmark_results.csv" \
  --out-raw "${CGAT_DIR}/generated/benchmark_results_raw.csv"

echo "[3/4] Update CGAT LaTeX tables from raw benchmark CSV..."
python3 "${CGAT_DIR}/tools/generate_cgat_tables.py" \
  --input "${CGAT_DIR}/generated/benchmark_results_raw.csv" \
  --outdir "${CGAT_DIR}/generated"

echo "[4/4] Build submission PDF..."
cd "${CGAT_DIR}"
pdflatex -interaction=nonstopmode -halt-on-error submission.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error submission.tex >/dev/null

echo "Done:"
echo "  - ${CGAT_DIR}/generated/benchmark_table.tex"
echo "  - ${CGAT_DIR}/generated/benchmark_full.tex"
echo "  - ${CGAT_DIR}/submission.pdf"

