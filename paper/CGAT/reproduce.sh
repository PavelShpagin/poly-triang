#!/usr/bin/env bash
set -euo pipefail

CGAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${CGAT_OUTDIR:-${CGAT_DIR}/generated}"
SKIP_PDF="${CGAT_SKIP_PDF:-0}"
KEEP_TEX_TEMP="${CGAT_KEEP_TEX_TEMP:-0}"

echo "[1/5] Build binaries..."
chmod +x "${CGAT_DIR}/build.sh" "${CGAT_DIR}/tools/"*.py
"${CGAT_DIR}/build.sh"

echo "[2/5] Sanity-check correctness (no crossings on sampled diagonals)..."
python3 "${CGAT_DIR}/tools/check_diagonal_validity.py" --n 500 --seed 0
python3 "${CGAT_DIR}/tools/check_diagonal_validity.py" --n 500 --seed 1

echo "[2b/5] Visual sanity-check: render + validate triangulations..."
FIG_DIR="${OUTDIR}/figures"
VIZ_N="${CGAT_VIZ_N:-200}"
VIZ_SEED="${CGAT_VIZ_SEED:-0}"

# Phase 3 (full triangulation) + validator (strong end-to-end correctness signal).
python3 "${CGAT_DIR}/tools/visualize_triangulation.py" \
  --bin "${CGAT_DIR}/bin/reflex_cli" \
  --mode triangulation \
  --outdir "${FIG_DIR}" \
  --n "${VIZ_N}" \
  --seed "${VIZ_SEED}"
python3 "${CGAT_DIR}/tools/validate_tri.py" --tri "${FIG_DIR}/triangulation_convex.tri"
python3 "${CGAT_DIR}/tools/validate_tri.py" --tri "${FIG_DIR}/triangulation_dent.tri"
python3 "${CGAT_DIR}/tools/validate_tri.py" --tri "${FIG_DIR}/triangulation_random.tri"

# Phase 2 (decomposition diagonals) for paper figures / debugging.
python3 "${CGAT_DIR}/tools/visualize_triangulation.py" \
  --mode decomposition \
  --diag-debug "${CGAT_DIR}/bin/diag_debug_cli" \
  --outdir "${FIG_DIR}" \
  --n "${VIZ_N}" \
  --seed "${VIZ_SEED}"

echo "[3/5] Run benchmarks (paper sizes)..."
TYPES="${CGAT_TYPES:-convex,dent,random,star}"
SIZES="${CGAT_SIZES:-500,1000,2000,5000,10000}"
RUNS="${CGAT_RUNS:-5}"
TIMEOUT="${CGAT_TIMEOUT:-10}"
PIN_CPU="${CGAT_PIN_CPU:-0}"

python3 "${CGAT_DIR}/tools/benchmark_cgat.py" \
  --types "${TYPES}" \
  --sizes "${SIZES}" \
  --runs "${RUNS}" \
  --timeout "${TIMEOUT}" \
  --pin-cpu "${PIN_CPU}" \
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
echo "  - ${FIG_DIR}/triangulation_*.svg"
echo "  - ${FIG_DIR}/triangulation_*.tri"
if [ "${SKIP_PDF}" != "1" ]; then
  echo "  - ${CGAT_DIR}/submission.pdf"
fi

