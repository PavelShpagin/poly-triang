#!/usr/bin/env python3
"""
Visualization script for polygon triangulation benchmarks.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (10, 6)

COLORS = {
    'earcut': '#e41a1c',           # Red - Optimized Earcut
    'earclip_naive': '#ff7f00',    # Orange - Naive O(n^2)
    'garey': '#377eb8',            # Blue - O(n log n)
    'hertel': '#4daf4a',           # Green - O(n + r log r)
}

LABELS = {
    'earcut': 'Earcut (optimized C++)',
    'earclip_naive': 'Ear Clipping naive O(n^2)',
    'garey': 'Garey O(n log n)',
    'hertel': 'Hertel-Mehlhorn O(n + r log r)',
}


def load_triangulation(filename):
    """Load triangulation from .tri file"""
    with open(filename, 'r') as f:
        # Skip comment lines
        line = f.readline().strip()
        while line.startswith('#'):
            line = f.readline().strip()
        n = int(line)
        
        vertices = []
        for _ in range(n):
            x, y = map(float, f.readline().strip().split())
            vertices.append((x, y))
        
        # Skip comment lines
        line = f.readline().strip()
        while line.startswith('#'):
            line = f.readline().strip()
        m = int(line)
        
        triangles = []
        for _ in range(m):
            v0, v1, v2 = map(int, f.readline().strip().split())
            triangles.append((v0, v1, v2))
    
    return np.array(vertices), triangles


def plot_triangulation(vertices, triangles, title, ax, color='#377eb8'):
    """Plot a single triangulation"""
    patches = []
    for tri in triangles:
        triangle = vertices[list(tri)]
        patches.append(MplPolygon(triangle, closed=True))
    
    p = PatchCollection(patches, alpha=0.4, facecolor=color, edgecolor='#333333', linewidth=0.5)
    ax.add_collection(p)
    
    poly_closed = np.vstack([vertices, vertices[0]])
    ax.plot(poly_closed[:, 0], poly_closed[:, 1], 'k-', linewidth=1.5)
    ax.scatter(vertices[:, 0], vertices[:, 1], c='black', s=20, zorder=5)
    
    ax.set_aspect('equal')
    ax.set_title(title)


def plot_triangulation_comparison(results_dir, output_dir):
    """Plot triangulations for all algorithms on sample polygons"""
    algorithms = list(COLORS.keys())
    poly_types = ['random_100', 'star_100']
    
    for poly_type in poly_types:
        fig, axes = plt.subplots(1, len(algorithms), figsize=(5 * len(algorithms), 5))
        if len(algorithms) == 1:
            axes = [axes]
        
        for idx, alg in enumerate(algorithms):
            tri_file = Path(results_dir) / f'{alg}_{poly_type}.tri'
            if not tri_file.exists():
                axes[idx].set_title(f'{alg} (not available)')
                continue
            
            vertices, triangles = load_triangulation(tri_file)
            plot_triangulation(vertices, triangles, 
                             f'{alg} ({len(triangles)} triangles)', 
                             axes[idx], COLORS[alg])
        
        plt.suptitle(f'Triangulation: {poly_type.replace("_", " ").title()}', fontsize=16)
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f'triangulation_{poly_type}.png', dpi=150, bbox_inches='tight')
        plt.close()


def plot_benchmark_times(df, output_dir):
    """Generate benchmark performance plots"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for alg in df['algorithm'].unique():
        alg_data = df[df['algorithm'] == alg].groupby('num_vertices')['time_ms'].mean().reset_index()
        label = LABELS.get(alg, alg)
        ax.plot(alg_data['num_vertices'], alg_data['time_ms'], 
               'o-', label=label, color=COLORS.get(alg, 'gray'), linewidth=2, markersize=6)
    
    ax.set_xlabel('Number of Vertices (N)')
    ax.set_ylabel('Time (ms)')
    ax.set_title('Triangulation Algorithm Performance Comparison')
    ax.legend(loc='upper left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'benchmark_times.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_by_polygon_type(df, output_dir):
    """Plot performance by polygon type"""
    poly_types = df['polygon'].str.extract(r'^([a-z]+)_')[0].unique()
    
    n_types = len(poly_types)
    fig, axes = plt.subplots(1, n_types, figsize=(5 * n_types, 5))
    if n_types == 1:
        axes = [axes]
    
    for idx, ptype in enumerate(poly_types):
        ax = axes[idx]
        subset = df[df['polygon'].str.startswith(ptype)]
        
        for alg in subset['algorithm'].unique():
            alg_data = subset[subset['algorithm'] == alg].sort_values('num_vertices')
            label = LABELS.get(alg, alg)
            ax.plot(alg_data['num_vertices'], alg_data['time_ms'], 
                   'o-', label=label, color=COLORS.get(alg, 'gray'), linewidth=2, markersize=5)
        
        ax.set_xlabel('N (vertices)')
        ax.set_ylabel('Time (ms)')
        ax.set_title(f'{ptype.title()} Polygons')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        if idx == n_types - 1:
            ax.legend(loc='upper left', fontsize=8)
    
    plt.suptitle('Algorithm Performance by Polygon Type', fontsize=14)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'benchmark_by_type.png', dpi=150, bbox_inches='tight')
    plt.close()


def generate_latex_table(df, output_path):
    """Generate LaTeX table for the paper"""
    pivot = df.pivot_table(
        values='time_ms', 
        index='polygon',
        columns='algorithm'
    ).round(4)
    
    algorithms = list(pivot.columns)
    
    with open(output_path, 'w') as f:
        f.write("% Auto-generated benchmark results table\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{Triangulation Algorithm Performance (time in ms)}\n")
        f.write("\\label{tab:benchmark_results}\n")
        f.write("\\begin{tabular}{l" + "r" * len(algorithms) + "}\n")
        f.write("\\toprule\n")
        f.write("Polygon & " + " & ".join(algorithms) + " \\\\\n")
        f.write("\\midrule\n")
        
        for polygon, row in pivot.iterrows():
            values = [f"{row.get(alg, '-'):.4f}" if pd.notna(row.get(alg)) else '-' 
                     for alg in algorithms]
            f.write(f"{polygon} & " + " & ".join(values) + " \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"LaTeX table saved to {output_path}")


def fit_scaling_law(n_values, time_values):
    """
    Fit scaling law T = a * n^b using log-log linear regression.
    Returns (a, b, r_squared).
    """
    log_n = np.log(n_values)
    log_t = np.log(time_values)
    
    # Linear regression in log-log space
    coeffs = np.polyfit(log_n, log_t, 1)
    b = coeffs[0]  # exponent
    a = np.exp(coeffs[1])  # coefficient
    
    # R-squared
    log_t_pred = coeffs[0] * log_n + coeffs[1]
    ss_res = np.sum((log_t - log_t_pred) ** 2)
    ss_tot = np.sum((log_t - np.mean(log_t)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    return a, b, r_squared


def plot_scaling_analysis(df, output_dir):
    """
    Analyze and plot scaling laws for each algorithm.
    Fits T = a * n^b and compares to theoretical complexity.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left plot: Empirical data with fitted curves
    ax1 = axes[0]
    
    # Right plot: Scaling exponents comparison
    ax2 = axes[1]
    
    scaling_results = {}
    
    for alg in df['algorithm'].unique():
        alg_data = df[df['algorithm'] == alg].groupby('num_vertices')['time_ms'].mean().reset_index()
        alg_data = alg_data.sort_values('num_vertices')
        
        n_vals = alg_data['num_vertices'].values
        t_vals = alg_data['time_ms'].values
        
        # Filter out any zeros or negatives
        mask = (n_vals > 0) & (t_vals > 0)
        n_vals = n_vals[mask]
        t_vals = t_vals[mask]
        
        if len(n_vals) < 3:
            continue
        
        # Fit scaling law
        a, b, r2 = fit_scaling_law(n_vals, t_vals)
        scaling_results[alg] = {'a': a, 'b': b, 'r2': r2}
        
        label = LABELS.get(alg, alg)
        color = COLORS.get(alg, 'gray')
        
        # Plot empirical data
        ax1.scatter(n_vals, t_vals, color=color, s=50, zorder=5)
        
        # Plot fitted curve
        n_fit = np.logspace(np.log10(n_vals.min()), np.log10(n_vals.max()) + 0.5, 100)
        t_fit = a * n_fit ** b
        ax1.plot(n_fit, t_fit, '-', color=color, linewidth=2, 
                label=f'{label}\n$T = {a:.2e} \\cdot n^{{{b:.2f}}}$ ($R^2={r2:.3f}$)')
    
    ax1.set_xlabel('Number of Vertices (N)')
    ax1.set_ylabel('Time (ms)')
    ax1.set_title('Empirical Scaling Laws')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Add reference lines for theoretical complexities
    n_ref = np.logspace(1, 5, 100)
    base_time = 0.001  # baseline for reference
    
    # O(n) reference
    ax1.plot(n_ref, base_time * n_ref / 10, '--', color='gray', alpha=0.5, linewidth=1, label='_O(n)')
    # O(n log n) reference  
    ax1.plot(n_ref, base_time * n_ref * np.log2(n_ref) / 100, '-.', color='gray', alpha=0.5, linewidth=1, label='_O(n log n)')
    # O(n^2) reference
    ax1.plot(n_ref, base_time * (n_ref ** 2) / 10000, ':', color='gray', alpha=0.5, linewidth=1, label='_O(n^2)')
    
    # Bar chart of scaling exponents
    algorithms = list(scaling_results.keys())
    exponents = [scaling_results[alg]['b'] for alg in algorithms]
    colors = [COLORS.get(alg, 'gray') for alg in algorithms]
    
    bars = ax2.bar(range(len(algorithms)), exponents, color=colors, alpha=0.7, edgecolor='black')
    ax2.set_xticks(range(len(algorithms)))
    ax2.set_xticklabels([LABELS.get(alg, alg).split()[0] for alg in algorithms], rotation=45, ha='right')
    ax2.set_ylabel('Scaling Exponent (b in T ~ n^b)')
    ax2.set_title('Empirical Scaling Exponents')
    ax2.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='O(n): b=1')
    ax2.axhline(y=2.0, color='red', linestyle='--', linewidth=2, label='O(n^2): b=2')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, exp in zip(bars, exponents):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05, 
                f'{exp:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'scaling_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Print scaling results
    print("\n  Scaling Analysis Results:")
    print("  " + "-" * 60)
    for alg, results in scaling_results.items():
        label = LABELS.get(alg, alg)
        print(f"  {label}:")
        print(f"    T = {results['a']:.4e} * n^{results['b']:.3f}")
        print(f"    R^2 = {results['r2']:.4f}")
    print("  " + "-" * 60)
    
    return scaling_results


def plot_extrapolation(df, output_dir, scaling_results):
    """
    Extrapolate performance to larger polygon sizes (millions of vertices).
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Extrapolation range: 10 to 10 million vertices
    n_extrap = np.logspace(1, 7, 200)
    
    for alg in scaling_results:
        a = scaling_results[alg]['a']
        b = scaling_results[alg]['b']
        
        # Extrapolated time
        t_extrap = a * n_extrap ** b
        
        label = LABELS.get(alg, alg)
        color = COLORS.get(alg, 'gray')
        
        ax.plot(n_extrap, t_extrap, '-', color=color, linewidth=2.5, label=label)
        
        # Plot actual data points
        alg_data = df[df['algorithm'] == alg].groupby('num_vertices')['time_ms'].mean().reset_index()
        ax.scatter(alg_data['num_vertices'], alg_data['time_ms'], color=color, s=80, 
                  edgecolor='black', linewidth=1, zorder=5)
    
    # Add reference lines
    ax.axhline(y=1000, color='gray', linestyle=':', alpha=0.7)
    ax.text(1.5e7, 1200, '1 second', fontsize=10, color='gray')
    ax.axhline(y=60000, color='gray', linestyle=':', alpha=0.7)
    ax.text(1.5e7, 72000, '1 minute', fontsize=10, color='gray')
    ax.axhline(y=3600000, color='gray', linestyle=':', alpha=0.7)
    ax.text(1.5e7, 4300000, '1 hour', fontsize=10, color='gray')
    
    # Add vertical lines for key sizes
    for n, label in [(1000, '1K'), (10000, '10K'), (100000, '100K'), (1000000, '1M'), (10000000, '10M')]:
        ax.axvline(x=n, color='lightgray', linestyle='-', alpha=0.3)
        ax.text(n, ax.get_ylim()[0] * 1.5, label, fontsize=9, ha='center', color='gray')
    
    ax.set_xlabel('Number of Vertices (N)', fontsize=12)
    ax.set_ylabel('Estimated Time (ms)', fontsize=12)
    ax.set_title('Extrapolated Performance: From Thousands to Millions of Vertices', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(10, 2e7)
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'extrapolation.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Print extrapolation estimates
    print("\n  Extrapolated Performance Estimates:")
    print("  " + "-" * 70)
    print(f"  {'Algorithm':<30} {'1K verts':<12} {'10K verts':<12} {'100K verts':<12} {'1M verts':<12}")
    print("  " + "-" * 70)
    for alg in scaling_results:
        a = scaling_results[alg]['a']
        b = scaling_results[alg]['b']
        label = LABELS.get(alg, alg).split('O(')[0].strip()
        
        t_1k = a * 1000 ** b
        t_10k = a * 10000 ** b
        t_100k = a * 100000 ** b
        t_1m = a * 1000000 ** b
        
        def format_time(ms):
            if ms < 1000:
                return f"{ms:.1f}ms"
            elif ms < 60000:
                return f"{ms/1000:.1f}s"
            elif ms < 3600000:
                return f"{ms/60000:.1f}min"
            else:
                return f"{ms/3600000:.1f}hr"
        
        print(f"  {label:<30} {format_time(t_1k):<12} {format_time(t_10k):<12} {format_time(t_100k):<12} {format_time(t_1m):<12}")
    print("  " + "-" * 70)


def generate_scaling_table(scaling_results, output_path):
    """Generate LaTeX table for scaling analysis."""
    with open(output_path, 'w') as f:
        f.write("% Auto-generated scaling analysis table\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{Empirical Scaling Law Analysis: $T = a \\cdot n^b$}\n")
        f.write("\\label{tab:scaling}\n")
        f.write("\\begin{tabular}{lcccc}\n")
        f.write("\\toprule\n")
        f.write("Algorithm & Theoretical & Empirical $b$ & Coefficient $a$ & $R^2$ \\\\\n")
        f.write("\\midrule\n")
        
        theoretical = {
            'earcut': 'Optimized C++',
            'earclip_naive': '$O(n^2)$, $b=2$',
            'garey': '$O(n \\log n)$, $b \\approx 1$',
            'hertel': '$O(n + r \\log r)$',
        }
        
        # Order: earclip_naive first to show true O(n^2), then others
        order = ['earclip_naive', 'earcut', 'garey', 'hertel']
        for alg in order:
            if alg not in scaling_results:
                continue
            results = scaling_results[alg]
            label = LABELS.get(alg, alg)
            if 'O(' in label:
                label = label.split('O(')[0].strip()
            theo = theoretical.get(alg, '-')
            f.write(f"{label} & {theo} & {results['b']:.3f} & {results['a']:.2e} & {results['r2']:.4f} \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"Scaling table saved to {output_path}")


def main():
    project_root = Path(__file__).parent.parent
    results_dir = project_root / 'results'
    output_dir = project_root / 'paper' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = results_dir / 'benchmark_results.csv'
    
    if not csv_path.exists():
        print(f"Error: Benchmark results not found at {csv_path}")
        sys.exit(1)
    
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} benchmark results")
    print(f"Polygon sizes: {sorted(df['num_vertices'].unique())}")
    
    print("\nGenerating visualizations...")
    
    print("  - Triangulation comparisons...")
    plot_triangulation_comparison(results_dir, output_dir)
    
    print("  - Benchmark performance plots...")
    plot_benchmark_times(df, output_dir)
    
    print("  - Performance by polygon type...")
    plot_by_polygon_type(df, output_dir)
    
    print("  - Scaling law analysis...")
    scaling_results = plot_scaling_analysis(df, output_dir)
    
    print("  - Extrapolation to large polygons...")
    plot_extrapolation(df, output_dir, scaling_results)
    
    print("  - LaTeX tables...")
    generate_latex_table(df, project_root / 'paper' / 'benchmark_table.tex')
    generate_scaling_table(scaling_results, project_root / 'paper' / 'scaling_table.tex')
    
    print(f"\nDone! Figures saved to {output_dir}")


if __name__ == '__main__':
    main()
