#!/usr/bin/env bash
# Quick paper benchmark - focuses on key comparisons
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
BIN_DIR="${BUILD_DIR}/bin"
RESULTS_DIR="${ROOT}/results"
POLY_DIR="${ROOT}/polygons/benchmark"
PAPER_GEN="${ROOT}/paper/generated"

# Smaller, faster configuration
REPEATS=3
SIZES="100 500 1000 2000 5000 10000 20000"
POLY_TYPES="convex random star spiral"
TIMEOUT_SEC=10

mkdir -p "${RESULTS_DIR}" "${POLY_DIR}" "${PAPER_GEN}"

echo "Quick Paper Benchmark"
echo "====================="

# Generate polygons
echo "[1/3] Generating polygons..."
python3 "${ROOT}/scripts/generate_polygons.py" --output "${POLY_DIR}" --sizes ${SIZES}

# Run benchmark
echo "[2/3] Running benchmarks..."
RAW_CSV="${RESULTS_DIR}/quick_benchmark_raw.csv"
echo "algorithm,polygon_type,num_vertices,num_reflex,time_ms,run" > "${RAW_CSV}"

run_algo() {
    local algo="$1"
    local cli="$2"
    local poly="$3"
    local ptype="$4"
    local n="$5"
    local run="$6"
    
    [ -x "${cli}" ] || return 0
    
    local output="${RESULTS_DIR}/tmp_${algo}.tri"
    if log="$(timeout ${TIMEOUT_SEC}s "${cli}" --input "${poly}" --output "${output}" 2>&1)"; then
        local time_ms reflex
        time_ms="$(echo "${log}" | grep -oP 'time_ms=\K[0-9.]+' || echo "0")"
        reflex="$(echo "${log}" | grep -oP 'reflex_count=\K[0-9]+' || echo "0")"
        echo "${algo},${ptype},${n},${reflex},${time_ms},${run}" >> "${RAW_CSV}"
    else
        # Timeout or failure: record as missing so tables show "--".
        echo "${algo},${ptype},${n},,,$run" >> "${RAW_CSV}"
    fi
    rm -f "${output}"
}

declare -A CLIS=(
    ["ours"]="${BIN_DIR}/reflex_cli"
    ["seidel"]="${BIN_DIR}/seidel_cli"
    ["garey"]="${BIN_DIR}/polypartition_mono_cli"
    ["hertel"]="${BIN_DIR}/polypartition_hm_cli"
)

for ptype in ${POLY_TYPES}; do
    for n in ${SIZES}; do
        poly="${POLY_DIR}/${ptype}_${n}.poly"
        [ -f "${poly}" ] || continue
        
        for run in $(seq 1 ${REPEATS}); do
            for algo in "${!CLIS[@]}"; do
                run_algo "${algo}" "${CLIS[$algo]}" "${poly}" "${ptype}" "${n}" "${run}"
            done
        done
        echo "  ${ptype} n=${n} done"
    done
done

# Generate tables
echo "[3/3] Generating LaTeX tables..."
python3 - "${RAW_CSV}" "${PAPER_GEN}" << 'PYLATEX'
import sys
import pandas as pd

raw_csv = sys.argv[1]
out_dir = sys.argv[2]

df = pd.read_csv(raw_csv)

# Aggregate
agg = df.groupby(['algorithm', 'polygon_type', 'num_vertices']).agg(
    time_mean=('time_ms', 'mean'),
    time_std=('time_ms', 'std'),
    reflex_mean=('num_reflex', 'mean'),
).reset_index()

NAMES = {
    'ours': r'\textbf{Ours}',
    'seidel': 'Seidel',
    'garey': 'Garey',
    'hertel': 'Hertel--Mehlhorn',
}
ORDER = ['ours', 'seidel', 'garey', 'hertel']

# Main table: random polygons
random_df = agg[agg['polygon_type'] == 'random'].copy()
if not random_df.empty:
    pivot = random_df.pivot_table(index='num_vertices', columns='algorithm', values='time_mean')
    # Reflex count is a property of the polygon; only our CLI reports it (baselines report 0).
    reflex_by_n = (
        random_df[random_df['algorithm'] == 'ours']
        .set_index('num_vertices')['reflex_mean']
    )
    cols = [c for c in ORDER if c in pivot.columns]
    pivot = pivot[cols]
    
    with open(f"{out_dir}/benchmark_table.tex", 'w') as f:
        f.write("\\begin{table}[t]\n\\centering\n")
        f.write("\\caption{Running time (ms) on random polygons.}\n")
        f.write("\\label{tab:benchmark}\n")
        f.write("\\begin{tabular}{rr" + "r" * len(cols) + "}\n")
        f.write("\\toprule\n")
        f.write("$n$ & $r$ & " + " & ".join(NAMES.get(c, c) for c in cols) + " \\\\\n")
        f.write("\\midrule\n")
        
        for n in sorted(pivot.index):
            row = pivot.loc[n]
            r = int(reflex_by_n.get(n, 0))
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
            f.write(f"{n:,} & {r} & " + " & ".join(vals) + " \\\\\n")
        
        f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

# Full appendix table
with open(f"{out_dir}/benchmark_full.tex", 'w') as f:
    f.write("\\begin{table*}[t]\n\\centering\n")
    f.write("\\caption{Running time (ms) across all polygon families.}\n")
    f.write("\\label{tab:benchmark-full}\n\\small\n")
    f.write("\\begin{tabular}{llr" + "r" * len(ORDER) + "}\n")
    f.write("\\toprule\n")
    f.write("Type & $n$ & $r$ & " + " & ".join(NAMES.get(c, c) for c in ORDER) + " \\\\\n")
    f.write("\\midrule\n")
    
    for ptype in ['convex', 'random', 'star', 'spiral']:
        ptype_df = agg[agg['polygon_type'] == ptype]
        if ptype_df.empty:
            continue
        
        pivot = ptype_df.pivot_table(index='num_vertices', columns='algorithm', values='time_mean')
        # Reflex count is a property of the polygon; only our CLI reports it (baselines report 0).
        reflex_by_n = (
            ptype_df[ptype_df['algorithm'] == 'ours']
            .set_index('num_vertices')['reflex_mean']
        )
        
        first_row = True
        for n in sorted(pivot.index):
            row = pivot.loc[n]
            r = int(reflex_by_n.get(n, 0))
            
            vals = []
            avail = [c for c in ORDER if c in row.index and not pd.isna(row.get(c))]
            min_val = row[avail].min() if avail else float('inf')
            
            for c in ORDER:
                if c not in pivot.columns or pd.isna(row.get(c)):
                    vals.append("--")
                elif row[c] == min_val:
                    vals.append(f"\\textbf{{{row[c]:.2f}}}")
                else:
                    vals.append(f"{row[c]:.2f}")
            
            type_col = ptype.capitalize() if first_row else ""
            f.write(f"{type_col} & {n:,} & {r} & " + " & ".join(vals) + " \\\\\n")
            first_row = False
        
        f.write("\\midrule\n")
    
    f.write("\\bottomrule\n\\end{tabular}\n\\end{table*}\n")

print("Tables generated!")
print(f"\nSummary (random polygons, mean time in ms):")
if not random_df.empty:
    summary = random_df.pivot_table(index='num_vertices', columns='algorithm', values='time_mean')
    print(summary.to_string())
PYLATEX

echo ""
echo "Done! Tables in ${PAPER_GEN}/"

