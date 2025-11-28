#!/usr/bin/env python3
"""
Rigorous visualization for polygon triangulation benchmark.
Creates two main figures:
1. Performance comparison across all polygon types
2. Triangulation visualization comparison
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['legend.fontsize'] = 10

COLORS = {
    'earcut': '#27ae60',      # Green
    'reflex': '#e74c3c',      # Red  
    'naive': '#7f8c8d',       # Gray
}

LABELS = {
    'earcut': 'Earcut (Mapbox C++)',
    'reflex': 'Reflex (Ours, C++)',
    'naive': 'Naive O(n^2) (Python)',
}

MARKERS = {
    'earcut': 'o',
    'reflex': 's',
    'naive': '^',
}


def fit_power_law(x, y):
    """Fit T = a * n^b using log-log linear regression."""
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


def load_triangulation(filename):
    """Load triangulation from .tri file."""
    with open(filename, 'r') as f:
        line = f.readline().strip()
        while line.startswith('#'):
            line = f.readline().strip()
        n = int(line)
        
        vertices = []
        for _ in range(n):
            x, y = map(float, f.readline().strip().split())
            vertices.append((x, y))
        
        line = f.readline().strip()
        while line.startswith('#'):
            line = f.readline().strip()
        m = int(line)
        
        triangles = []
        for _ in range(m):
            parts = f.readline().strip().split()
            v0, v1, v2 = int(parts[0]), int(parts[1]), int(parts[2])
            triangles.append((v0, v1, v2))
    
    return np.array(vertices), triangles


def figure1_performance_comparison(df, output_dir):
    """
    Figure 1: Performance comparison across all polygon types.
    Two subplots: (a) by polygon type, (b) overall with error bars.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    # (a) Performance by polygon type
    ax1 = axes[0]
    
    poly_types = ['convex', 'random', 'star', 'spiral']
    type_colors = {'convex': '#3498db', 'random': '#e74c3c', 'star': '#f39c12', 'spiral': '#9b59b6'}
    
    for ptype in poly_types:
        pdata = df[(df['polygon_type'] == ptype) & (df['algorithm'] == 'earcut')].sort_values('num_vertices')
        if len(pdata) > 0:
            ax1.plot(pdata['num_vertices'], pdata['time_ms'], 'o--', 
                    color=type_colors[ptype], alpha=0.6, linewidth=1.5, markersize=5,
                    label=f'Earcut - {ptype}')
        
        pdata = df[(df['polygon_type'] == ptype) & (df['algorithm'] == 'reflex')].sort_values('num_vertices')
        if len(pdata) > 0:
            ax1.plot(pdata['num_vertices'], pdata['time_ms'], 's-', 
                    color=type_colors[ptype], linewidth=2, markersize=6,
                    label=f'Reflex - {ptype}')
    
    ax1.set_xlabel('Number of Vertices (n)')
    ax1.set_ylabel('Time (ms)')
    ax1.set_title('(a) Performance by Polygon Type')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend(loc='upper left', fontsize=8, ncol=2)
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_xlim(50, 200000)
    
    # (b) Average performance with complexity fits
    ax2 = axes[1]
    
    print("\n  Complexity Analysis (T = a * n^b):")
    print("  " + "-" * 50)
    
    for alg in ['reflex', 'earcut', 'naive']:
        if alg not in df['algorithm'].unique():
            continue
        
        alg_df = df[df['algorithm'] == alg]
        grouped = alg_df.groupby('num_vertices')['time_ms'].agg(['mean', 'std']).reset_index()
        grouped = grouped.sort_values('num_vertices')
        
        x = grouped['num_vertices'].values
        y = grouped['mean'].values
        yerr = grouped['std'].values
        
        ax2.errorbar(x, y, yerr=yerr, fmt=MARKERS[alg] + '-', 
                    label=LABELS[alg], color=COLORS[alg], 
                    linewidth=2, markersize=7, capsize=3)
        
        # Fit power law
        a, b, r2 = fit_power_law(x, y)
        if a is not None:
            x_fit = np.logspace(np.log10(x.min()), np.log10(x.max()), 50)
            ax2.plot(x_fit, a * (x_fit ** b), '--', color=COLORS[alg], 
                    alpha=0.4, linewidth=1.5)
            print(f"  {LABELS[alg]:30s}: b = {b:.3f}, R^2 = {r2:.4f}")
    
    print("  " + "-" * 50)
    
    # Reference lines
    x_ref = np.logspace(1.5, 5.5, 100)
    scale = 1e-5
    ax2.plot(x_ref, scale * x_ref, 'k:', alpha=0.3, linewidth=1.5, label=r'$O(n)$')
    ax2.plot(x_ref, scale * x_ref * np.log2(x_ref), 'k-.', alpha=0.3, linewidth=1.5, label=r'$O(n \log n)$')
    ax2.plot(x_ref, scale * (x_ref ** 2) / 100, 'k--', alpha=0.3, linewidth=1.5, label=r'$O(n^2)$')
    
    ax2.set_xlabel('Number of Vertices (n)')
    ax2.set_ylabel('Time (ms)')
    ax2.set_title('(b) Average Performance (All Types)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend(loc='upper left', fontsize=9)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_xlim(50, 200000)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'figure1_performance.png', dpi=200, bbox_inches='tight')
    plt.savefig(output_dir / 'figure1_performance.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"\n  Saved figure1_performance.png")


def figure2_triangulation_comparison(results_dir, output_dir):
    """
    Figure 2: Visual comparison of triangulations on sample polygons.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    
    test_cases = [
        ('random_500', 'Random Polygon (500 vertices)'),
        ('star_500', 'Star Polygon (500 vertices)'),
    ]
    
    algorithms = ['earcut', 'reflex']
    
    for row, (poly_name, title) in enumerate(test_cases):
        for col, alg in enumerate(algorithms):
            ax = axes[row, col]
            tri_file = results_dir / f'{alg}_{poly_name}.tri'
            
            if not tri_file.exists():
                ax.text(0.5, 0.5, f'{alg}\nNot available', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12)
                ax.axis('off')
                continue
            
            vertices, triangles = load_triangulation(tri_file)
            
            # Draw triangles
            patches = []
            for tri in triangles:
                triangle = vertices[list(tri)]
                patches.append(MplPolygon(triangle, closed=True))
            
            p = PatchCollection(patches, alpha=0.5, facecolor=COLORS[alg], 
                              edgecolor='#2c3e50', linewidth=0.4)
            ax.add_collection(p)
            
            # Draw polygon boundary
            poly_closed = np.vstack([vertices, vertices[0]])
            ax.plot(poly_closed[:, 0], poly_closed[:, 1], 'k-', linewidth=1.5)
            
            ax.set_aspect('equal')
            subtitle = f'{LABELS[alg]}\n{len(triangles)} triangles'
            ax.set_title(subtitle, fontsize=11)
            ax.axis('off')
    
    # Add row labels
    fig.text(0.02, 0.75, 'Random\nPolygon', ha='center', va='center', fontsize=12, fontweight='bold', rotation=90)
    fig.text(0.02, 0.25, 'Star\nPolygon', ha='center', va='center', fontsize=12, fontweight='bold', rotation=90)
    
    plt.suptitle('Triangulation Comparison (500 vertices)', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0.03, 0, 1, 0.96])
    plt.savefig(output_dir / 'figure2_triangulation.png', dpi=200, bbox_inches='tight')
    plt.savefig(output_dir / 'figure2_triangulation.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  Saved figure2_triangulation.png")


def generate_latex_table(df, output_path):
    """Generate LaTeX table with complete results - no NaN."""
    
    # Pivot table
    pivot = df.groupby(['num_vertices', 'algorithm'])['time_ms'].mean().unstack()
    
    with open(output_path, 'w') as f:
        f.write("% Auto-generated benchmark results table\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{Average Triangulation Time (ms) by Polygon Size}\n")
        f.write("\\label{tab:performance}\n")
        f.write("\\begin{tabular}{r")
        
        cols = ['reflex', 'earcut', 'naive']
        cols = [c for c in cols if c in pivot.columns]
        f.write("c" * len(cols))
        f.write("c}\n")  # Extra column for speedup
        f.write("\\toprule\n")
        
        headers = [c.title() for c in cols] + ['Speedup']
        f.write("$n$ & " + " & ".join(headers) + " \\\\\n")
        f.write("\\midrule\n")
        
        for size in sorted(pivot.index):
            row = [f"{size:,}"]
            earcut_val = pivot.loc[size, 'earcut'] if 'earcut' in pivot.columns else np.nan
            reflex_val = pivot.loc[size, 'reflex'] if 'reflex' in pivot.columns else np.nan
            
            for col in cols:
                val = pivot.loc[size, col] if col in pivot.columns else np.nan
                if pd.isna(val):
                    row.append("--")
                elif val < 0.1:
                    row.append(f"{val:.4f}")
                elif val < 1:
                    row.append(f"{val:.3f}")
                elif val < 100:
                    row.append(f"{val:.2f}")
                elif val < 10000:
                    row.append(f"{val:.1f}")
                else:
                    row.append(f"{val:.0f}")
            
            # Speedup
            if not pd.isna(earcut_val) and not pd.isna(reflex_val) and reflex_val > 0:
                speedup = earcut_val / reflex_val
                if speedup >= 1:
                    row.append(f"\\textbf{{{speedup:.1f}x}}")
                else:
                    row.append(f"{speedup:.2f}x")
            else:
                row.append("--")
            
            f.write(" & ".join(row) + " \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"  LaTeX table saved to {output_path}")


def analyze_complexity(df, output_path):
    """Analyze and report complexity for each algorithm."""
    
    results = []
    
    for alg in df['algorithm'].unique():
        alg_df = df[df['algorithm'] == alg]
        grouped = alg_df.groupby('num_vertices')['time_ms'].mean().reset_index()
        grouped = grouped.sort_values('num_vertices')
        
        x = grouped['num_vertices'].values
        y = grouped['time_ms'].values
        
        a, b, r2 = fit_power_law(x, y)
        if a is not None:
            results.append({
                'Algorithm': LABELS.get(alg, alg),
                'Empirical b': b,
                'R^2': r2,
                'Coefficient a': a
            })
    
    with open(output_path, 'w') as f:
        f.write("Complexity Analysis: T = a * n^b\n")
        f.write("=" * 60 + "\n\n")
        
        for r in results:
            f.write(f"{r['Algorithm']}:\n")
            f.write(f"  T = {r['Coefficient a']:.2e} * n^{r['Empirical b']:.3f}\n")
            f.write(f"  R^2 = {r['R^2']:.4f}\n")
            
            b = r['Empirical b']
            if b < 1.15:
                interp = "Near-linear O(n)"
            elif b < 1.4:
                interp = "Between O(n) and O(n log n)"
            elif b < 1.7:
                interp = "Near O(n log n) to O(n sqrt(n))"
            else:
                interp = "Near quadratic O(n^2)"
            f.write(f"  Interpretation: {interp}\n\n")
        
        f.write("\nNOTE: The 'Reflex' algorithm shows high variance across polygon types.\n")
        f.write("- Convex polygons: O(n) via fan triangulation\n")
        f.write("- Non-convex: O(n * k) where k depends on spatial distribution\n")
        f.write("- Worst case: O(n^2) for pathological inputs\n")
    
    print(f"  Complexity analysis saved to {output_path}")


def main():
    project_root = Path(__file__).parent.parent
    results_dir = project_root / 'results'
    output_dir = project_root / 'publication' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = results_dir / 'rigorous_benchmark.csv'
    
    if not csv_path.exists():
        print(f"Error: Results not found at {csv_path}")
        sys.exit(1)
    
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} benchmark results")
    print(f"Algorithms: {df['algorithm'].unique().tolist()}")
    print(f"Sizes: {sorted(df['num_vertices'].unique())}")
    print(f"Polygon types: {df['polygon_type'].unique().tolist()}")
    
    print("\nGenerating figures...")
    
    figure1_performance_comparison(df, output_dir)
    figure2_triangulation_comparison(results_dir, output_dir)
    
    print("\n  Generating tables...")
    generate_latex_table(df, project_root / 'publication' / 'results_table.tex')
    analyze_complexity(df, project_root / 'publication' / 'complexity_analysis.txt')
    
    print(f"\nDone! Figures saved to {output_dir}")


if __name__ == '__main__':
    main()
