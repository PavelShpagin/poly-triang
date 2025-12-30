#!/usr/bin/env python3
"""
Visualization script for polygon triangulation methods comparison.
Focuses on comparing our reflex algorithm against baselines.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (12, 7)

COLORS = {
    'earcut': '#e41a1c',           # Red
    'reflex': '#984ea3',           # Purple - our method
    'earclip_naive': '#ff7f00',    # Orange
    'garey': '#377eb8',            # Blue
    'hertel': '#4daf4a',           # Green
    'kirkpatrick_seidel': '#a65628',  # Brown
}

LABELS = {
    'earcut': 'Earcut (optimized C++)',
    'reflex': 'Reflex O(n + r log r) [Ours]',
    'earclip_naive': 'Ear Clipping naive O(n^2)',
    'garey': 'Garey O(n log n)',
    'hertel': 'Hertel-Mehlhorn O(n + r log r)',
    'kirkpatrick_seidel': 'Kirkpatrick-Seidel (randomized) O(n log* n)',
}

MARKERS = {
    'earcut': 'o',
    'reflex': 's',
    'earclip_naive': '^',
    'garey': 'D',
    'hertel': 'v',
    'kirkpatrick_seidel': 'P',
}


def power_law(x, a, b):
    return a * (x ** b)


def log_power_law(x, log_a, b):
    return log_a + b * np.log(x)


def plot_benchmark_comparison(df, output_dir):
    """Main comparison plot: all algorithms on log-log scale."""
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Order: our method first, then others
    order = ['reflex', 'kirkpatrick_seidel', 'earcut', 'earclip_naive', 'garey', 'hertel']
    
    for alg in order:
        if alg not in df['algorithm'].unique():
            continue
        alg_data = df[df['algorithm'] == alg].groupby('num_vertices')['time_ms'].mean().reset_index()
        alg_data = alg_data.sort_values('num_vertices')
        
        label = LABELS.get(alg, alg)
        ax.plot(alg_data['num_vertices'], alg_data['time_ms'], 
               marker=MARKERS.get(alg, 'o'), label=label, 
               color=COLORS.get(alg, 'gray'), linewidth=2.5, markersize=8,
               markeredgecolor='white', markeredgewidth=1)
    
    ax.set_xlabel('Number of Vertices (N)')
    ax.set_ylabel('Time (ms)')
    ax.set_title('Triangulation Algorithm Performance Comparison')
    ax.legend(loc='upper left', framealpha=0.9)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'methods_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_scaling_analysis(df, output_dir):
    """Scaling law analysis with power law fits."""
    fig, ax = plt.subplots(figsize=(12, 7))
    
    scaling_results = {}
    order = ['reflex', 'kirkpatrick_seidel', 'earcut', 'earclip_naive', 'garey', 'hertel']
    
    for alg in order:
        if alg not in df['algorithm'].unique():
            continue
        alg_data = df[df['algorithm'] == alg].groupby('num_vertices')['time_ms'].mean().reset_index()
        alg_data = alg_data[alg_data['time_ms'] > 0].sort_values('num_vertices')
        
        if len(alg_data) < 3:
            continue
        
        x_data = alg_data['num_vertices'].values
        y_data = alg_data['time_ms'].values
        
        try:
            popt, _ = curve_fit(log_power_law, x_data, np.log(y_data), p0=[np.log(y_data[0]), 1])
            log_a_fit, b_fit = popt
            a_fit = np.exp(log_a_fit)
            
            y_pred = power_law(x_data, a_fit, b_fit)
            ss_res = np.sum((np.log(y_data) - np.log(y_pred))**2)
            ss_tot = np.sum((np.log(y_data) - np.mean(np.log(y_data)))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            scaling_results[alg] = {'a': a_fit, 'b': b_fit, 'r2': r_squared}
            
            # Plot empirical data
            label = f'{LABELS.get(alg, alg).split("[")[0].strip()}'
            ax.plot(x_data, y_data, marker=MARKERS.get(alg, 'o'), 
                   label=f'{label} (b={b_fit:.2f})', 
                   color=COLORS.get(alg, 'gray'), linewidth=2, markersize=7,
                   markeredgecolor='white', markeredgewidth=0.5)
            
            # Plot fit line
            x_fit = np.logspace(np.log10(x_data.min()), np.log10(x_data.max()), 50)
            ax.plot(x_fit, power_law(x_fit, a_fit, b_fit), '--', 
                   color=COLORS.get(alg, 'gray'), alpha=0.5, linewidth=1.5)
            
        except Exception as e:
            print(f"  Warning: Could not fit {alg}: {e}")
    
    # Add theoretical reference lines
    n_ref = np.logspace(1, 4.5, 100)
    ref_scale = 0.001
    ax.plot(n_ref, ref_scale * n_ref**2, 'k:', alpha=0.4, linewidth=1.5, label=r'$O(n^2)$ reference')
    ax.plot(n_ref, ref_scale * n_ref * np.log2(n_ref), 'k-.', alpha=0.4, linewidth=1.5, label=r'$O(n \log n)$ reference')
    ax.plot(n_ref, ref_scale * n_ref, 'k--', alpha=0.4, linewidth=1.5, label=r'$O(n)$ reference')
    
    ax.set_xlabel('Number of Vertices (N)')
    ax.set_ylabel('Time (ms)')
    ax.set_title('Scaling Law Analysis: Empirical Complexity')
    ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'methods_scaling.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return scaling_results


def plot_by_polygon_type(df, output_dir):
    """Performance breakdown by polygon type."""
    poly_types = df['polygon'].str.extract(r'^([a-z]+)_')[0].unique()
    
    n_types = len(poly_types)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    order = ['reflex', 'earcut', 'earclip_naive', 'garey', 'hertel']
    
    for idx, ptype in enumerate(poly_types[:4]):
        ax = axes[idx]
        subset = df[df['polygon'].str.startswith(ptype)]
        
        for alg in order:
            if alg not in subset['algorithm'].unique():
                continue
            alg_data = subset[subset['algorithm'] == alg].sort_values('num_vertices')
            label = LABELS.get(alg, alg).split('O(')[0].strip()
            ax.plot(alg_data['num_vertices'], alg_data['time_ms'], 
                   marker=MARKERS.get(alg, 'o'), label=label, 
                   color=COLORS.get(alg, 'gray'), linewidth=2, markersize=6)
        
        ax.set_xlabel('N (vertices)')
        ax.set_ylabel('Time (ms)')
        ax.set_title(f'{ptype.title()} Polygons')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        if idx == 1:
            ax.legend(loc='upper left', fontsize=8)
    
    plt.suptitle('Algorithm Performance by Polygon Type', fontsize=14)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'methods_by_type.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_reflex_impact(df, output_dir):
    """Show how reflex count affects our algorithm's performance."""
    reflex_data = df[df['algorithm'] == 'reflex'].copy()
    
    if reflex_data.empty or 'num_reflex' not in reflex_data.columns:
        print("  Skipping reflex impact plot (no data)")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: Time vs reflex count
    ax1.scatter(reflex_data['num_reflex'], reflex_data['time_ms'], 
               c=reflex_data['num_vertices'], cmap='viridis', s=80, alpha=0.7)
    ax1.set_xlabel('Number of Reflex Vertices (r)')
    ax1.set_ylabel('Time (ms)')
    ax1.set_title('Reflex Algorithm: Time vs Reflex Count')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis', 
                               norm=plt.Normalize(vmin=reflex_data['num_vertices'].min(),
                                                  vmax=reflex_data['num_vertices'].max()))
    cbar = plt.colorbar(sm, ax=ax1)
    cbar.set_label('Total Vertices (n)')
    
    # Right: Reflex ratio by polygon type
    reflex_data['reflex_ratio'] = reflex_data['num_reflex'] / reflex_data['num_vertices']
    reflex_data['poly_type'] = reflex_data['polygon'].str.extract(r'^([a-z]+)_')[0]
    
    for ptype in reflex_data['poly_type'].unique():
        pdata = reflex_data[reflex_data['poly_type'] == ptype].sort_values('num_vertices')
        ax2.plot(pdata['num_vertices'], pdata['reflex_ratio'], 'o-', label=ptype.title(), linewidth=2)
    
    ax2.set_xlabel('Number of Vertices (N)')
    ax2.set_ylabel('Reflex Ratio (r/n)')
    ax2.set_title('Reflex Vertex Ratio by Polygon Type')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'reflex_impact.png', dpi=150, bbox_inches='tight')
    plt.close()


def generate_latex_table(df, scaling_results, output_path):
    """Generate LaTeX table for the publication."""
    with open(output_path, 'w') as f:
        f.write("% Auto-generated methods comparison table\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{Empirical Scaling Analysis of Triangulation Algorithms}\n")
        f.write("\\label{tab:methods_scaling}\n")
        f.write("\\begin{tabular}{lccccc}\n")
        f.write("\\toprule\n")
        f.write("Algorithm & Theoretical & Empirical $b$ & $R^2$ & 1K (ms) & 10K (ms) \\\\\n")
        f.write("\\midrule\n")
        
        theoretical = {
            'reflex': '$O(n + r \\log r)$',
            'earcut': 'Optimized',
            'earclip_naive': '$O(n^2)$',
            'garey': '$O(n \\log n)$',
            'hertel': '$O(n + r \\log r)$',
        }
        
        order = ['reflex', 'earcut', 'earclip_naive', 'garey', 'hertel']
        
        for alg in order:
            if alg not in scaling_results:
                continue
            
            results = scaling_results[alg]
            name = LABELS.get(alg, alg).split('O(')[0].strip()
            if '[Ours]' in LABELS.get(alg, ''):
                name = '\\textbf{' + name.replace('[Ours]', '').strip() + '} (Ours)'
            
            theo = theoretical.get(alg, '-')
            
            # Get actual times for 1K and 10K
            alg_data = df[df['algorithm'] == alg]
            t_1k = alg_data[alg_data['num_vertices'] == 1000]['time_ms'].mean()
            t_10k = alg_data[alg_data['num_vertices'] == 10000]['time_ms'].mean()
            
            t_1k_str = f"{t_1k:.3f}" if not np.isnan(t_1k) else "-"
            t_10k_str = f"{t_10k:.2f}" if not np.isnan(t_10k) else "-"
            
            f.write(f"{name} & {theo} & {results['b']:.3f} & {results['r2']:.3f} & {t_1k_str} & {t_10k_str} \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"  LaTeX table saved to {output_path}")


def main():
    project_root = Path(__file__).parent.parent
    results_dir = project_root / 'results'
    output_dir = project_root / 'publication' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = results_dir / 'methods_benchmark.csv'
    
    if not csv_path.exists():
        print(f"Error: Benchmark results not found at {csv_path}")
        print("Run scripts/benchmark_methods.sh first")
        sys.exit(1)
    
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} benchmark results")
    print(f"Algorithms: {df['algorithm'].unique().tolist()}")
    print(f"Polygon sizes: {sorted(df['num_vertices'].unique())}")
    
    print("\nGenerating visualizations...")
    
    print("  - Main comparison plot...")
    plot_benchmark_comparison(df, output_dir)
    
    print("  - Scaling analysis...")
    scaling_results = plot_scaling_analysis(df, output_dir)
    
    print("  - Performance by polygon type...")
    plot_by_polygon_type(df, output_dir)
    
    print("  - Reflex impact analysis...")
    plot_reflex_impact(df, output_dir)
    
    print("  - LaTeX table...")
    generate_latex_table(df, scaling_results, project_root / 'publication' / 'methods_table.tex')
    
    # Print scaling summary
    print("\n  Scaling Analysis Summary:")
    print("  " + "-" * 60)
    for alg, res in scaling_results.items():
        name = LABELS.get(alg, alg).split('O(')[0].strip()
        print(f"  {name:30s} b={res['b']:.3f}  R^2={res['r2']:.4f}")
    print("  " + "-" * 60)
    
    print(f"\nDone! Figures saved to {output_dir}")


if __name__ == '__main__':
    main()

