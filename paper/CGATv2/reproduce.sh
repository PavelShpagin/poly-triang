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
OUTDIR="${CGAT_OUTDIR:-${CGAT_DIR}/generated}"
SKIP_PDF="${CGAT_SKIP_PDF:-0}"
KEEP_TEX_TEMP="${CGAT_KEEP_TEX_TEMP:-0}"

python3 "${CGAT_DIR}/tools/benchmark_cgat.py" \
  --types convex,dent,random,star \
  --sizes 500,1000,2000,5000,10000 \
  --runs 5 \
  --timeout 10 \
  --bin-dir "${CGAT_DIR}/bin" \
  --out-csv "${OUTDIR}/benchmark_results.csv" \
  --out-raw "${OUTDIR}/benchmark_results_raw.csv"

echo "[4/5] Update CGAT LaTeX tables from raw benchmark CSV..."
python3 "${CGAT_DIR}/tools/generate_cgat_tables.py" \
  --input "${OUTDIR}/benchmark_results_raw.csv" \
  --outdir "${OUTDIR}"

echo "[5/5] Build submission PDF..."
cd "${CGAT_DIR}"
if [ "${SKIP_PDF}" != "1" ]; then
  pdflatex -interaction=nonstopmode -halt-on-error submission.tex
  pdflatex -interaction=nonstopmode -halt-on-error submission.tex
fi

if [ "${KEEP_TEX_TEMP}" != "1" ]; then
  echo "Cleaning LaTeX temporary files..."
  rm -f submission.aux submission.log submission.out submission.spl submission.toc \
        submission.bbl submission.blg submission.fls submission.fdb_latexmk \
        submission.lof submission.lot submission.synctex.gz 2>/dev/null || true
fi

echo "Done:"
echo "  - ${OUTDIR}/benchmark_table.tex"
echo "  - ${OUTDIR}/benchmark_full.tex"
if [ "${SKIP_PDF}" != "1" ]; then
  echo "  - ${CGAT_DIR}/submission.pdf"
fi

