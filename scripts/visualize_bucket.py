#!/usr/bin/env python3
"""
Visualization for Bucket Triangulation benchmark.
Demonstrates O(n) amortized complexity.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['legend.fontsize'] = 11

COLORS = {
    'bucket': '#e74c3c',      # Red - our method
    'earcut': '#27ae60',      # Green - baseline
}

LABELS = {
    'bucket': 'Bucket (Ours) - Amortized O(n)',
    'earcut': 'Earcut (Mapbox)',
}


def fit_power_law(x, y):
    """Fit T = a * n^b."""
    try:
        valid = (y > 0) & (x > 0)
        if np.sum(valid) < 3:
            return None, None, None
        log_x = np.log(x[valid])
        log_y = np.log(y[valid])
        coeffs = np.polyfit(log_x, log_y, 1)
        b = coeffs[0]
        a = np.exp(coeffs[1])
        y_pred = a * (x[valid] ** b)
        ss_res = np.sum((y[valid] - y_pred)**2)
        ss_tot = np.sum((y[valid] - np.mean(y[valid]))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return a, b, r2
    except:
        return None, None, None


def main_figure(df, output_dir):
    """
    Main figure: Performance comparison showing O(n) scaling.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    # (a) Log-log performance plot
    ax1 = axes[0]
    
    scaling_results = {}
    
    for alg in ['bucket', 'earcut']:
        if alg not in df['algorithm'].unique():
            continue
        
        alg_df = df[df['algorithm'] == alg]
        grouped = alg_df.groupby('num_vertices')['time_ms'].agg(['mean', 'std']).reset_index()
        grouped = grouped.sort_values('num_vertices')
        
        x = grouped['num_vertices'].values
        y = grouped['mean'].values
        yerr = grouped['std'].values
        
        ax1.errorbar(x, y, yerr=yerr, fmt='o-', label=LABELS[alg], 
                    color=COLORS[alg], linewidth=2.5, markersize=8, capsize=4)
        
        # Fit power law
        a, b, r2 = fit_power_law(x, y)
        if a is not None:
            scaling_results[alg] = {'a': a, 'b': b, 'r2': r2}
            x_fit = np.logspace(np.log10(x.min()), np.log10(x.max()), 50)
            ax1.plot(x_fit, a * (x_fit ** b), '--', color=COLORS[alg], 
                    alpha=0.5, linewidth=2)
    
    # Reference lines
    x_ref = np.logspace(3, 6.5, 100)
    scale = 1e-5
    ax1.plot(x_ref, scale * x_ref, 'k:', alpha=0.4, linewidth=2, label=r'$O(n)$')
    ax1.plot(x_ref, scale * x_ref * np.log2(x_ref) / 5, 'k-.', alpha=0.4, linewidth=2, label=r'$O(n \log n)$')
    ax1.plot(x_ref, scale * (x_ref ** 1.5) / 50, 'k--', alpha=0.4, linewidth=2, label=r'$O(n^{1.5})$')
    
    ax1.set_xlabel('Number of Vertices (n)')
    ax1.set_ylabel('Time (ms)')
    ax1.set_title('(a) Performance Scaling to 1M Vertices')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_xlim(500, 2000000)
    
    # (b) Time per vertex (should be constant for O(n))
    ax2 = axes[1]
    
    for alg in ['bucket', 'earcut']:
        if alg not in df['algorithm'].unique():
            continue
        
        alg_df = df[df['algorithm'] == alg]
        grouped = alg_df.groupby('num_vertices')['time_ms'].mean().reset_index()
        grouped = grouped.sort_values('num_vertices')
        
        x = grouped['num_vertices'].values
        y = grouped['time_ms'].values
        
        # Time per vertex in microseconds
        time_per_vertex = (y / x) * 1000  # Convert to microseconds
        
        ax2.plot(x, time_per_vertex, 'o-', label=LABELS[alg], 
                color=COLORS[alg], linewidth=2.5, markersize=8)
    
    ax2.set_xlabel('Number of Vertices (n)')
    ax2.set_ylabel('Time per Vertex (microseconds)')
    ax2.set_title('(b) Amortized Cost per Vertex')
    ax2.set_xscale('log')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_xlim(500, 2000000)
    
    # Add horizontal reference line for O(n) behavior
    ax2.axhline(y=ax2.get_ylim()[0] + (ax2.get_ylim()[1] - ax2.get_ylim()[0]) * 0.3, 
               color='gray', linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'figure1_scaling.png', dpi=200, bbox_inches='tight')
    plt.savefig(output_dir / 'figure1_scaling.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  Saved figure1_scaling.png")
    
    return scaling_results


def speedup_figure(df, output_dir):
    """
    Figure showing speedup by polygon type.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Compute speedup
    bucket_df = df[df['algorithm'] == 'bucket'].copy()
    earcut_df = df[df['algorithm'] == 'earcut'].copy()
    
    merged = pd.merge(bucket_df, earcut_df, on=['polygon_type', 'num_vertices'], 
                     suffixes=('_bucket', '_earcut'))
    merged['speedup'] = merged['time_ms_earcut'] / merged['time_ms_bucket']
    
    poly_types = merged['polygon_type'].unique()
    colors = {'convex': '#3498db', 'random': '#e74c3c', 'star': '#f39c12', 'spiral': '#9b59b6'}
    
    for ptype in poly_types:
        pdata = merged[merged['polygon_type'] == ptype].sort_values('num_vertices')
        if len(pdata) > 0:
            ax.plot(pdata['num_vertices'], pdata['speedup'], 'o-', 
                   label=ptype.title(), color=colors.get(ptype, 'gray'), 
                   linewidth=2, markersize=7)
    
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, label='Equal performance')
    ax.set_xlabel('Number of Vertices (n)')
    ax.set_ylabel('Speedup (Earcut time / Bucket time)')
    ax.set_title('Speedup of Bucket Algorithm over Earcut')
    ax.set_xscale('log')
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(500, 2000000)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'figure2_speedup.png', dpi=200, bbox_inches='tight')
    plt.savefig(output_dir / 'figure2_speedup.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  Saved figure2_speedup.png")


def generate_latex_table(df, scaling_results, output_path):
    """Generate LaTeX table with results."""
    
    pivot = df.groupby(['num_vertices', 'algorithm'])['time_ms'].mean().unstack()
    
    with open(output_path, 'w') as f:
        f.write("% Auto-generated benchmark results\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{Triangulation Performance (ms) - Bucket vs Earcut}\n")
        f.write("\\label{tab:bucket_performance}\n")
        f.write("\\begin{tabular}{rccc}\n")
        f.write("\\toprule\n")
        f.write("$n$ & Bucket (Ours) & Earcut & Speedup \\\\\n")
        f.write("\\midrule\n")
        
        for size in sorted(pivot.index):
            bucket_val = pivot.loc[size, 'bucket'] if 'bucket' in pivot.columns else np.nan
            earcut_val = pivot.loc[size, 'earcut'] if 'earcut' in pivot.columns else np.nan
            
            def fmt(v):
                if pd.isna(v): return "--"
                if v < 1: return f"{v:.3f}"
                if v < 100: return f"{v:.2f}"
                if v < 10000: return f"{v:.1f}"
                return f"{v:.0f}"
            
            speedup = earcut_val / bucket_val if not pd.isna(bucket_val) and bucket_val > 0 else np.nan
            speedup_str = f"{speedup:.2f}x" if not pd.isna(speedup) else "--"
            if speedup and speedup > 1:
                speedup_str = f"\\textbf{{{speedup_str}}}"
            
            f.write(f"{size:,} & {fmt(bucket_val)} & {fmt(earcut_val)} & {speedup_str} \\\\\n")
        
        f.write("\\midrule\n")
        
        # Add scaling exponents
        if 'bucket' in scaling_results:
            b = scaling_results['bucket']['b']
            f.write(f"Scaling $b$ & {b:.3f} & ")
        if 'earcut' in scaling_results:
            b = scaling_results['earcut']['b']
            f.write(f"{b:.3f} & -- \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"  LaTeX table saved to {output_path}")


def main():
    project_root = Path(__file__).parent.parent
    results_dir = project_root / 'results'
    output_dir = project_root / 'experiments' / 'bucket' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = results_dir / 'bucket_benchmark.csv'
    
    if not csv_path.exists():
        print(f"Error: Results not found at {csv_path}")
        sys.exit(1)
    
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} benchmark results")
    print(f"Algorithms: {df['algorithm'].unique().tolist()}")
    print(f"Max vertices: {df['num_vertices'].max():,}")
    
    print("\nGenerating figures...")
    
    scaling_results = main_figure(df, output_dir)
    speedup_figure(df, output_dir)
    
    print("\n  Scaling Analysis:")
    print("  " + "-" * 50)
    for alg, res in scaling_results.items():
        print(f"  {LABELS.get(alg, alg):35s}: b = {res['b']:.3f}, R^2 = {res['r2']:.4f}")
    print("  " + "-" * 50)
    
    generate_latex_table(df, scaling_results, project_root / 'experiments' / 'bucket' / 'results_table.tex')
    
    print(f"\nDone! Figures saved to {output_dir}")


if __name__ == '__main__':
    main()

