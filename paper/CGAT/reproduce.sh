#!/usr/bin/env bash
set -euo pipefail

CGTA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Preferred variables are CGTA_*. Legacy CGAT_* variables are still accepted.
OUTDIR="${CGTA_OUTDIR:-${CGAT_OUTDIR:-${CGTA_DIR}/generated}}"
SKIP_PDF="${CGTA_SKIP_PDF:-${CGAT_SKIP_PDF:-0}}"
KEEP_TEX_TEMP="${CGTA_KEEP_TEX_TEMP:-${CGAT_KEEP_TEX_TEMP:-0}}"

echo "[1/5] Build binaries..."
chmod +x "${CGTA_DIR}/build.sh" "${CGTA_DIR}/tools/"*.py
"${CGTA_DIR}/build.sh"

echo "[2/5] Sanity-check correctness (no crossings on sampled diagonals)..."
python3 "${CGTA_DIR}/tools/check_diagonal_validity.py" --n 500 --seed 0
python3 "${CGTA_DIR}/tools/check_diagonal_validity.py" --n 500 --seed 1

echo "[2b/5] Visual sanity-check: render + validate triangulations..."
FIG_DIR="${OUTDIR}/figures"
VIZ_N="${CGTA_VIZ_N:-${CGAT_VIZ_N:-200}}"
VIZ_SEED="${CGTA_VIZ_SEED:-${CGAT_VIZ_SEED:-0}}"

# Phase 3 (full triangulation) + validator (strong end-to-end correctness signal).
python3 "${CGTA_DIR}/tools/visualize_triangulation.py" \
  --bin "${CGTA_DIR}/bin/reflex_cli" \
  --mode triangulation \
  --outdir "${FIG_DIR}" \
  --n "${VIZ_N}" \
  --seed "${VIZ_SEED}"
python3 "${CGTA_DIR}/tools/validate_tri.py" --tri "${FIG_DIR}/triangulation_convex.tri"
python3 "${CGTA_DIR}/tools/validate_tri.py" --tri "${FIG_DIR}/triangulation_dent.tri"
python3 "${CGTA_DIR}/tools/validate_tri.py" --tri "${FIG_DIR}/triangulation_random.tri"

# Phase 2 (decomposition diagonals) for paper figures / debugging.
python3 "${CGTA_DIR}/tools/visualize_triangulation.py" \
  --mode decomposition \
  --diag-debug "${CGTA_DIR}/bin/diag_debug_cli" \
  --outdir "${FIG_DIR}" \
  --n "${VIZ_N}" \
  --seed "${VIZ_SEED}"

echo "[3/5] Run benchmarks (paper sizes)..."
TYPES="${CGTA_TYPES:-${CGAT_TYPES:-convex,dent,random,star}}"
SIZES="${CGTA_SIZES:-${CGAT_SIZES:-500,1000,2000,5000,10000}}"
RUNS="${CGTA_RUNS:-${CGAT_RUNS:-5}}"
TIMEOUT="${CGTA_TIMEOUT:-${CGAT_TIMEOUT:-10}}"
PIN_CPU="${CGTA_PIN_CPU:-${CGAT_PIN_CPU:-0}}"

python3 "${CGTA_DIR}/tools/benchmark_cgat.py" \
  --types "${TYPES}" \
  --sizes "${SIZES}" \
  --runs "${RUNS}" \
  --timeout "${TIMEOUT}" \
  --pin-cpu "${PIN_CPU}" \
  --bin-dir "${CGTA_DIR}/bin" \
  --out-csv "${OUTDIR}/benchmark_results.csv" \
  --out-raw "${OUTDIR}/benchmark_results_raw.csv"

echo "[4/5] Update CGTA LaTeX tables from raw benchmark CSV..."
python3 "${CGTA_DIR}/tools/generate_cgat_tables.py" \
  --input "${OUTDIR}/benchmark_results_raw.csv" \
  --outdir "${OUTDIR}"

echo "[5/5] Build submission PDF..."
cd "${CGTA_DIR}"
if [ "${SKIP_PDF}" != "1" ]; then
  pdflatex -interaction=nonstopmode -halt-on-error submission.tex
  pdflatex -interaction=nonstopmode -halt-on-error submission.tex
  if [ -f cover_letter.tex ]; then
    pdflatex -interaction=nonstopmode -halt-on-error cover_letter.tex
    pdflatex -interaction=nonstopmode -halt-on-error cover_letter.tex
  fi
fi

if [ "${KEEP_TEX_TEMP}" != "1" ]; then
  echo "Cleaning LaTeX temporary files..."
  rm -f submission.aux submission.log submission.out submission.spl submission.toc \
        submission.bbl submission.blg submission.fls submission.fdb_latexmk \
        submission.lof submission.lot submission.synctex.gz \
        cover_letter.aux cover_letter.log cover_letter.out cover_letter.toc \
        cover_letter.fls cover_letter.fdb_latexmk cover_letter.synctex.gz 2>/dev/null || true
fi

echo "Done:"
echo "  - ${OUTDIR}/benchmark_table.tex"
echo "  - ${OUTDIR}/benchmark_full.tex"
echo "  - ${FIG_DIR}/triangulation_*.svg"
echo "  - ${FIG_DIR}/triangulation_*.tri"
if [ "${SKIP_PDF}" != "1" ]; then
  echo "  - ${CGTA_DIR}/submission.pdf"
  if [ -f "${CGTA_DIR}/cover_letter.pdf" ]; then
    echo "  - ${CGTA_DIR}/cover_letter.pdf"
  fi
fi

