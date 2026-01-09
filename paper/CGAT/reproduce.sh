#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CGAT_DIR="${ROOT}/paper/CGAT"

echo "[1/4] Build binaries (Release)..."
mkdir -p "${ROOT}/build"
cd "${ROOT}/build"
cmake -DCMAKE_BUILD_TYPE=Release .. >/dev/null
make -j4 reflex_cli polypartition_mono_cli polypartition_hm_cli seidel_cli >/dev/null

echo "[2/4] Run benchmarks (paper sizes)..."
cd "${ROOT}"
python3 scripts/benchmark_paper.py \
  --types convex,dent,random,star \
  --sizes 100,500,1000,2000,5000,10000 \
  --runs 5 \
  --timeout 10 \
  --out-csv "${CGAT_DIR}/generated/benchmark_results.csv" \
  --out-raw-csv "${CGAT_DIR}/generated/benchmark_results_raw.csv"

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

