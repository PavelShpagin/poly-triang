#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${ROOT}/bin"

mkdir -p "${BIN_DIR}"

CXX="${CXX:-g++}"
CC="${CC:-gcc}"

CXXFLAGS=(-O3 -DNDEBUG -std=c++17)
CFLAGS=(-O3 -DNDEBUG -std=gnu99 -D_DEFAULT_SOURCE)

SEIDEL_DIR="${ROOT}/code/baselines/Seidel"
SEIDEL_DEFS=(-DSEGSIZE=20000)
SEIDEL_INC=(-I "${SEIDEL_DIR}/Seidel")

echo "[0/3] Build Seidel objects (for linking into ours)..."
"${CC}" "${CFLAGS[@]}" "${SEIDEL_DEFS[@]}" "${SEIDEL_INC[@]}" \
  -c "${SEIDEL_DIR}/Seidel/construct.c" -o "${BIN_DIR}/seidel_construct.o"
"${CC}" "${CFLAGS[@]}" "${SEIDEL_DEFS[@]}" "${SEIDEL_INC[@]}" \
  -c "${SEIDEL_DIR}/Seidel/misc.c" -o "${BIN_DIR}/seidel_misc.o"
"${CC}" "${CFLAGS[@]}" "${SEIDEL_DEFS[@]}" "${SEIDEL_INC[@]}" \
  -c "${SEIDEL_DIR}/Seidel/monotone.c" -o "${BIN_DIR}/seidel_monotone.o"
"${CC}" "${CFLAGS[@]}" "${SEIDEL_DEFS[@]}" "${SEIDEL_INC[@]}" \
  -c "${SEIDEL_DIR}/Seidel/tri.c" -o "${BIN_DIR}/seidel_tri.o"

echo "[1/3] Build ours..."
"${CXX}" "${CXXFLAGS[@]}" \
  "${SEIDEL_DEFS[@]}" \
  -I "${ROOT}/code/ours" \
  "${SEIDEL_INC[@]}" \
  "${ROOT}/code/ours/reflex_cli.cpp" \
  "${BIN_DIR}/seidel_construct.o" \
  "${BIN_DIR}/seidel_misc.o" \
  "${BIN_DIR}/seidel_monotone.o" \
  "${BIN_DIR}/seidel_tri.o" \
  -lm \
  -o "${BIN_DIR}/reflex_cli"

echo "[1b/3] Build diag_debug_cli (for correctness checks)..."
"${CXX}" "${CXXFLAGS[@]}" \
  "${SEIDEL_DEFS[@]}" \
  -I "${ROOT}/code/ours" \
  "${SEIDEL_INC[@]}" \
  "${ROOT}/code/ours/diag_debug_cli.cpp" \
  "${BIN_DIR}/seidel_construct.o" \
  "${BIN_DIR}/seidel_misc.o" \
  "${BIN_DIR}/seidel_monotone.o" \
  "${BIN_DIR}/seidel_tri.o" \
  -lm \
  -o "${BIN_DIR}/diag_debug_cli"

echo "[2/3] Build PolyPartition baselines..."
"${CXX}" "${CXXFLAGS[@]}" \
  -I "${ROOT}/code/baselines/polypartition/src" \
  "${ROOT}/code/cli/polypartition_mono_cli.cpp" \
  "${ROOT}/code/baselines/polypartition/src/polypartition.cpp" \
  -o "${BIN_DIR}/polypartition_mono_cli"

"${CXX}" "${CXXFLAGS[@]}" \
  "${ROOT}/code/cli/cgal_hm_cli.cpp" \
  -lgmp -lmpfr \
  -o "${BIN_DIR}/cgal_hm_cli"

echo "[3/3] Build Seidel baseline..."
"${CC}" "${CFLAGS[@]}" \
  "${SEIDEL_DEFS[@]}" \
  "${SEIDEL_INC[@]}" \
  "${ROOT}/code/cli/seidel_cli.c" \
  "${SEIDEL_DIR}/Seidel/construct.c" \
  "${SEIDEL_DIR}/Seidel/misc.c" \
  "${SEIDEL_DIR}/Seidel/monotone.c" \
  "${SEIDEL_DIR}/Seidel/tri.c" \
  -lm \
  -o "${BIN_DIR}/seidel_cli"

echo "Built:"
echo "  - ${BIN_DIR}/reflex_cli"
echo "  - ${BIN_DIR}/diag_debug_cli"
echo "  - ${BIN_DIR}/polypartition_mono_cli"
echo "  - ${BIN_DIR}/cgal_hm_cli"
echo "  - ${BIN_DIR}/seidel_cli"

