#!/usr/bin/env python3
"""Analyze C++ benchmark results and generate paper tables."""
import pandas as pd
from pathlib import Path

# Load results
csv_path = Path(__file__).parent.parent / 'results' / 'cpp_benchmark.csv'
df = pd.read_csv(csv_path)

# Add reflex ratio where available
df['r_ratio'] = df['r'] / df['n']

# Group by polygon type, n, and algorithm, compute mean
summary = df.groupby(['polygon_type', 'n', 'algorithm']).agg({
    'time_ms': 'mean',
    'r': 'first',
    'r_ratio': 'first'
}).reset_index()

print("=" * 80)
print("COMPREHENSIVE BENCHMARK RESULTS")
print("=" * 80)

# Create pivot table for each polygon type
for ptype in ['convex', 'spiral', 'random', 'star']:
    pdata = summary[summary['polygon_type'] == ptype]
    if pdata.empty:
        continue
    
    pivot = pdata.pivot_table(values='time_ms', index='n', columns='algorithm')
    
    # Get r values
    r_vals = pdata[pdata['algorithm'] == 'reflex'].set_index('n')['r']
    
    print(f"\n{ptype.upper()} POLYGONS:")
    print("-" * 60)
    
    for n in sorted(pivot.index):
        reflex_t = pivot.loc[n, 'reflex'] if 'reflex' in pivot.columns else float('nan')
        earcut_t = pivot.loc[n, 'earcut'] if 'earcut' in pivot.columns else float('nan')
        polytri_t = pivot.loc[n, 'polytri'] if 'polytri' in pivot.columns else float('nan')
        r = r_vals.get(n, 0)
        
        speedup_earcut = earcut_t / reflex_t if reflex_t > 0 else 0
        speedup_polytri = polytri_t / reflex_t if reflex_t > 0 else 0
        
        print(f"  n={n:>6}  r={r:>5}  Ours={reflex_t:>8.3f}ms  "
              f"Earcut={earcut_t:>8.3f}ms  Speedup={speedup_earcut:>6.2f}x")

# Generate LaTeX table
print("\n" + "=" * 80)
print("LATEX TABLE FOR PAPER")
print("=" * 80)

latex_output = []
latex_output.append("\\begin{table}[t]")
latex_output.append("\\centering")
latex_output.append("\\caption{Running times (ms) comparing our $O(n + r \\log r)$ algorithm against Earcut.}")
latex_output.append("\\label{tab:benchmark}")
latex_output.append("\\begin{tabular}{llrrrrr}")
latex_output.append("\\toprule")
latex_output.append("\\textbf{Type} & $n$ & $r$ & \\textbf{Ours} & \\textbf{Earcut} & \\textbf{Speedup} \\\\")
latex_output.append("\\midrule")

for ptype in ['convex', 'spiral', 'star']:
    pdata = summary[summary['polygon_type'] == ptype]
    pivot = pdata.pivot_table(values='time_ms', index='n', columns='algorithm')
    r_vals = pdata[pdata['algorithm'] == 'reflex'].set_index('n')['r']
    
    for n in [1000, 5000, 10000, 50000]:
        if n not in pivot.index:
            continue
        reflex_t = pivot.loc[n, 'reflex']
        earcut_t = pivot.loc[n, 'earcut']
        r = int(r_vals.get(n, 0))
        speedup = earcut_t / reflex_t
        
        ptype_display = ptype.title()
        latex_output.append(f"{ptype_display} & {n:,} & {r:,} & {reflex_t:.2f} & {earcut_t:.2f} & {speedup:.1f}$\\times$ \\\\")

latex_output.append("\\bottomrule")
latex_output.append("\\end{tabular}")
latex_output.append("\\end{table}")

print("\n".join(latex_output))

# Save LaTeX to file
latex_path = Path(__file__).parent / 'benchmark_table.tex'
with open(latex_path, 'w') as f:
    f.write("\n".join(latex_output))
print(f"\nLaTeX saved to {latex_path}")

