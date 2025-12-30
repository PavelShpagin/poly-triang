#!/usr/bin/env bash
# ============================================================================
# PAPER BENCHMARK: One-shot script for reproducible experimental results
# ============================================================================
# This script:
#   1. Builds all triangulation CLIs (ours + baselines)
#   2. Generates polygon families (convex, random, star, spiral, CGAL random)
#   3. Runs all methods with multiple repeats
#   4. Outputs raw CSV and aggregated CSV
#   5. Generates LaTeX tables for the paper
#
# Usage: ./scripts/paper_benchmark.sh
# Output:
#   - results/paper_benchmark_raw.csv     (all runs)
#   - results/paper_benchmark.csv         (aggregated means/stdev)
#   - paper/generated/benchmark_table.tex (main paper table)
#   - paper/generated/benchmark_full.tex  (appendix table)
# ============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/benchmark"
PAPER_GEN="${ROOT}/paper/generated"
SCRIPTS_DIR="${ROOT}/scripts"

# Configuration
REPEATS=5
# Paper evaluation uses up to 20k for reproducible end-to-end runs.
SIZES="100 500 1000 2000 5000 10000 20000"
POLY_TYPES="convex random star spiral"
TIMEOUT_SEC=10

mkdir -p "${RESULTS_DIR}" "${POLY_DIR}" "${PAPER_GEN}"

echo "============================================================================"
echo "PAPER BENCHMARK: Polygon Triangulation Performance Evaluation"
echo "============================================================================"
echo ""

# ----------------------------------------------------------------------------
# Step 1: Build all executables
# ----------------------------------------------------------------------------
echo "[1/5] Building all triangulation CLIs..."
cmake -S "${ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 || {
    echo "CMake configuration failed. Check CMakeLists.txt"
    exit 1
}

TARGETS="reflex_cli earcut_cli polytri_cli polypartition_mono_cli polypartition_hm_cli seidel_cli"
for target in ${TARGETS}; do
    cmake --build "${BUILD_DIR}" --target "${target}" -j4 2>/dev/null || {
        echo "  Warning: Failed to build ${target} (may be optional)"
    }
done

# Check for CGAL-based generator (optional)
HAVE_CGAL=false
if cmake --build "${BUILD_DIR}" --target cgal_gen_cli -j4 2>/dev/null; then
    HAVE_CGAL=true
    echo "  CGAL generator available"
fi

echo "  Build complete"
echo ""

# ----------------------------------------------------------------------------
# Step 2: Generate polygon families
# ----------------------------------------------------------------------------
echo "[2/5] Generating polygon families..."
python3 "${SCRIPTS_DIR}/generate_polygons.py" --output "${POLY_DIR}" --sizes ${SIZES}

# Generate CGAL random_simple polygons if available
if [ "${HAVE_CGAL}" = true ]; then
    echo "  Generating CGAL random_simple polygons..."
    for n in ${SIZES}; do
        for seed in 1 2 3; do
            out="${POLY_DIR}/cgal_random_${n}_s${seed}.poly"
            if [ -f "${out}" ]; then
                continue
            fi
            "${BIN_DIR}/cgal_gen_cli" --type random_simple --n "${n}" --seed "${seed}" \
                --output "${out}" 2>/dev/null || true
        done
    done
fi
echo "  Polygon generation complete"
echo ""

# ----------------------------------------------------------------------------
# Step 3: Run benchmarks
# ----------------------------------------------------------------------------
echo "[3/5] Running benchmarks (${REPEATS} repeats per configuration)..."

RAW_CSV="${RESULTS_DIR}/paper_benchmark_raw.csv"
echo "algorithm,polygon_type,num_vertices,num_reflex,time_ms,run" > "${RAW_CSV}"

run_algo() {
    local algo="$1"
    local cli="$2"
    local poly="$3"
    local ptype="$4"
    local n="$5"
    local run="$6"
    
    if [ ! -x "${cli}" ]; then
        return
    fi
    
    local output="${RESULTS_DIR}/tmp_${algo}.tri"
    if log="$(timeout ${TIMEOUT_SEC}s "${cli}" --input "${poly}" --output "${output}" 2>&1)"; then
        local time_ms
        time_ms="$(echo "${log}" | grep -oP 'time_ms=\K[0-9.]+' || echo "0")"
        local reflex
        reflex="$(echo "${log}" | grep -oP 'reflex_count=\K[0-9]+' || echo "0")"
        echo "${algo},${ptype},${n},${reflex},${time_ms},${run}" >> "${RAW_CSV}"
    fi
    rm -f "${output}"
}

# Define CLI paths
declare -A CLIS=(
    ["reflex"]="${BIN_DIR}/reflex_cli"
    ["seidel"]="${BIN_DIR}/seidel_cli"
    ["garey"]="${BIN_DIR}/polypartition_mono_cli"
    ["hertel"]="${BIN_DIR}/polypartition_hm_cli"
    ["earcut"]="${BIN_DIR}/earcut_cli"
    ["polytri"]="${BIN_DIR}/polytri_cli"
)

total_runs=0
for ptype in ${POLY_TYPES}; do
    for n in ${SIZES}; do
        poly="${POLY_DIR}/${ptype}_${n}.poly"
        [ -f "${poly}" ] || continue
        
        for run in $(seq 1 ${REPEATS}); do
            for algo in "${!CLIS[@]}"; do
                run_algo "${algo}" "${CLIS[$algo]}" "${poly}" "${ptype}" "${n}" "${run}"
                ((total_runs++))
            done
        done
        printf "\r  Processed: %s n=%d (%d runs)     " "${ptype}" "${n}" "${total_runs}"
    done
done

# Run CGAL random if available
if [ "${HAVE_CGAL}" = true ]; then
    for n in ${SIZES}; do
        for seed in 1 2 3; do
            poly="${POLY_DIR}/cgal_random_${n}_s${seed}.poly"
            [ -f "${poly}" ] || continue
            
            for run in $(seq 1 ${REPEATS}); do
                for algo in "${!CLIS[@]}"; do
                    run_algo "${algo}" "${CLIS[$algo]}" "${poly}" "cgal_random" "${n}" "${run}"
                    ((total_runs++))
                done
            done
        done
        printf "\r  Processed: cgal_random n=%d (%d runs)     " "${n}" "${total_runs}"
    done
fi

echo ""
echo "  Benchmark runs complete: ${total_runs} total"
echo ""

# ----------------------------------------------------------------------------
# Step 4: Aggregate results
# ----------------------------------------------------------------------------
echo "[4/5] Aggregating results..."

AGG_CSV="${RESULTS_DIR}/paper_benchmark.csv"
python3 - "${RAW_CSV}" "${AGG_CSV}" << 'PYAGG'
import sys
import pandas as pd

raw_csv = sys.argv[1]
agg_csv = sys.argv[2]

df = pd.read_csv(raw_csv)

# Aggregate by (algorithm, polygon_type, num_vertices)
agg = df.groupby(['algorithm', 'polygon_type', 'num_vertices']).agg(
    time_mean=('time_ms', 'mean'),
    time_std=('time_ms', 'std'),
    time_min=('time_ms', 'min'),
    time_max=('time_ms', 'max'),
    reflex_mean=('num_reflex', 'mean'),
    runs=('run', 'count')
).reset_index()

agg.to_csv(agg_csv, index=False)
print(f"  Aggregated {len(df)} raw measurements into {len(agg)} configurations")
PYAGG

echo ""

# ----------------------------------------------------------------------------
# Step 5: Generate LaTeX tables
# ----------------------------------------------------------------------------
echo "[5/5] Generating LaTeX tables..."

python3 - "${AGG_CSV}" "${PAPER_GEN}" << 'PYLATEX'
import sys
import pandas as pd
import numpy as np

agg_csv = sys.argv[1]
out_dir = sys.argv[2]

df = pd.read_csv(agg_csv)

# Algorithm display names and order
NAMES = {
    'reflex': r'\textbf{Ours}',
    'seidel': 'Seidel',
    'garey': 'Garey',
    'hertel': 'Hertel--Mehlhorn',
    'earcut': 'Earcut',
    'polytri': 'PolyTri',
}
ORDER = ['reflex', 'seidel', 'garey', 'hertel', 'earcut', 'polytri']

# --- Main table: CGAL random_simple polygons only ---
random_df = df[df['polygon_type'] == 'cgal_random'].copy()
if not random_df.empty:
    pivot = random_df.pivot_table(
        index='num_vertices',
        columns='algorithm',
        values='time_mean',
        aggfunc='mean'
    )
    
    # Reorder columns
    cols = [c for c in ORDER if c in pivot.columns]
    pivot = pivot[cols]
    
    with open(f"{out_dir}/benchmark_table.tex", 'w') as f:
        f.write("% Auto-generated from paper_benchmark.csv\n")
        f.write("\\begin{table}[t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Running time (ms) on CGAL random\\_simple polygons.}\n")
        f.write("\\label{tab:benchmark}\n")
        f.write("\\begin{tabular}{r" + "r" * len(cols) + "}\n")
        f.write("\\toprule\n")
        f.write("$n$ & " + " & ".join(NAMES.get(c, c) for c in cols) + " \\\\\n")
        f.write("\\midrule\n")
        
        for n in sorted(pivot.index):
            row = pivot.loc[n]
            vals = []
            min_val = row.min()
            for c in cols:
                v = row[c]
                if pd.isna(v):
                    vals.append("--")
                elif v == min_val:
                    vals.append(f"\\textbf{{{v:.2f}}}")
                else:
                    vals.append(f"{v:.2f}")
            f.write(f"{n:,} & " + " & ".join(vals) + " \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    print(f"  Generated {out_dir}/benchmark_table.tex")

# --- Full appendix table: all polygon types ---
with open(f"{out_dir}/benchmark_full.tex", 'w') as f:
    f.write("% Auto-generated from paper_benchmark.csv\n")
    f.write("\\begin{table*}[t]\n")
    f.write("\\centering\n")
    f.write("\\caption{Running time (ms) across all polygon families.}\n")
    f.write("\\label{tab:benchmark-full}\n")
    f.write("\\small\n")
    f.write("\\begin{tabular}{llr" + "r" * len(ORDER) + "}\n")
    f.write("\\toprule\n")
    f.write("Type & $n$ & $r$ & " + " & ".join(NAMES.get(c, c) for c in ORDER) + " \\\\\n")
    f.write("\\midrule\n")
    
    for ptype in ['convex', 'cgal_random', 'star', 'spiral']:
        ptype_df = df[df['polygon_type'] == ptype]
        if ptype_df.empty:
            continue
        
        pivot = ptype_df.pivot_table(
            index='num_vertices',
            columns='algorithm',
            values='time_mean',
            aggfunc='mean'
        )
        reflex_by_n = ptype_df.groupby('num_vertices')['reflex_mean'].first()
        
        first_row = True
        for n in sorted(pivot.index):
            row = pivot.loc[n]
            r = int(reflex_by_n.get(n, 0))
            
            vals = []
            avail = [c for c in ORDER if c in row.index and not pd.isna(row[c])]
            min_val = row[avail].min() if avail else float('inf')
            
            for c in ORDER:
                if c not in pivot.columns or pd.isna(row.get(c)):
                    vals.append("--")
                elif row[c] == min_val:
                    vals.append(f"\\textbf{{{row[c]:.2f}}}")
                else:
                    vals.append(f"{row[c]:.2f}")
            
        type_col = ("Random" if ptype == "cgal_random" else ptype.replace('_', ' ').title()) if first_row else ""
            f.write(f"{type_col} & {n:,} & {r} & " + " & ".join(vals) + " \\\\\n")
            first_row = False
        
        f.write("\\midrule\n")
    
    f.write("\\bottomrule\n")
    f.write("\\end{tabular}\n")
    f.write("\\end{table*}\n")

print(f"  Generated {out_dir}/benchmark_full.tex")

# --- Summary statistics ---
print("\n  Summary by algorithm (all polygon types):")
summary = df.groupby('algorithm')['time_mean'].agg(['mean', 'std']).round(3)
for alg in ORDER:
    if alg in summary.index:
        print(f"    {NAMES.get(alg, alg):20s}: mean={summary.loc[alg, 'mean']:.3f}ms, std={summary.loc[alg, 'std']:.3f}ms")
PYLATEX

echo ""
echo "============================================================================"
echo "BENCHMARK COMPLETE"
echo "============================================================================"
echo "  Raw data:      ${RAW_CSV}"
echo "  Aggregated:    ${AGG_CSV}"
echo "  LaTeX table:   ${PAPER_GEN}/benchmark_table.tex"
echo "  LaTeX full:    ${PAPER_GEN}/benchmark_full.tex"
echo ""
echo "To include in paper, add to triangulation_note.tex:"
echo "  \\input{generated/benchmark_table}"
echo "============================================================================"

