#!/usr/bin/env python3
"""
Comprehensive benchmark comparing:
- Our O(n + r log r) algorithm
- Garey et al. O(n log n) baseline
"""
import sys
import time
import csv
import math
import random
from pathlib import Path

from triangulate_correct import triangulate as triangulate_ours, gen_convex, gen_star, gen_random_star, gen_spiral
from garey_baseline import triangulate_garey


def benchmark_algo(algo, pts, warmup=1, runs=3):
    """Benchmark a triangulation algorithm."""
    for _ in range(warmup):
        algo(pts)
    
    times = []
    r = 0
    for _ in range(runs):
        t0 = time.perf_counter()
        _, r = algo(pts)
        times.append((time.perf_counter() - t0) * 1000)
    
    return sum(times) / len(times), r


def run_comparison(gen, ptype, sizes, polygons_per_size=20, timing_runs=3):
    """Run comparison for a polygon type."""
    results = []
    
    for n in sizes:
        ours_times = []
        garey_times = []
        r_vals = []
        
        for seed in range(polygons_per_size):
            pts = gen(n, seed)
            
            t_ours, r = benchmark_algo(triangulate_ours, pts, warmup=1, runs=timing_runs)
            t_garey, _ = benchmark_algo(triangulate_garey, pts, warmup=1, runs=timing_runs)
            
            ours_times.append(t_ours)
            garey_times.append(t_garey)
            r_vals.append(r)
        
        avg_ours = sum(ours_times) / len(ours_times)
        avg_garey = sum(garey_times) / len(garey_times)
        avg_r = sum(r_vals) / len(r_vals)
        
        speedup = avg_garey / avg_ours if avg_ours > 0 else float('inf')
        
        results.append({
            'type': ptype,
            'n': n,
            'r': int(avg_r),
            'r_ratio': avg_r / n,
            'ours_ms': avg_ours,
            'garey_ms': avg_garey,
            'speedup': speedup,
            'samples': polygons_per_size,
        })
        
        print(f"  {ptype:<12} n={n:>6} r={int(avg_r):>5} ({avg_r/n:>5.1%}) "
              f"Ours={avg_ours:>8.2f}ms  Garey={avg_garey:>8.2f}ms  "
              f"Speedup={speedup:>6.2f}x", file=sys.stderr)
    
    return results


def main():
    print("=" * 85, file=sys.stderr)
    print("POLYGON TRIANGULATION BENCHMARK - O(n + r log r) vs O(n log n)", file=sys.stderr)
    print("=" * 85, file=sys.stderr)
    
    # Configuration
    sizes = [100, 500, 1000, 2000, 5000, 10000]
    polygons_per_size = 20  # Statistical stability
    timing_runs = 3
    
    polygon_types = [
        ('convex', gen_convex),
        ('spiral', gen_spiral),
        ('random_star', gen_random_star),
        ('star', gen_star),
    ]
    
    all_results = []
    
    for ptype, gen in polygon_types:
        print(f"\n{ptype.upper()}:", file=sys.stderr)
        results = run_comparison(gen, ptype, sizes, polygons_per_size, timing_runs)
        all_results.extend(results)
    
    # Save to CSV
    output_path = Path(__file__).parent / 'benchmark_comparison.csv'
    with open(output_path, 'w', newline='') as f:
        fieldnames = ['type', 'n', 'r', 'r_ratio', 'ours_ms', 'garey_ms', 'speedup', 'samples']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    
    print(f"\nResults saved to {output_path}", file=sys.stderr)
    
    # Print summary tables
    print("\n" + "=" * 85, file=sys.stderr)
    print("SUMMARY: Speedup over Garey et al.", file=sys.stderr)
    print("=" * 85, file=sys.stderr)
    
    for ptype in ['convex', 'spiral', 'random_star', 'star']:
        type_results = [r for r in all_results if r['type'] == ptype]
        if type_results:
            r10k = next((r for r in type_results if r['n'] == 10000), None)
            if r10k:
                print(f"  {ptype:<12}: {r10k['speedup']:.1f}x speedup at n=10000 (r={r10k['r']})", file=sys.stderr)
    
    # Generate LaTeX table for paper
    latex_path = Path(__file__).parent / 'benchmark_table.tex'
    with open(latex_path, 'w') as f:
        f.write("% Auto-generated benchmark table\n")
        f.write("\\begin{table}[t]\n")
        f.write("\\centering\n")
        f.write("\\caption{Running times (ms) comparing our $O(n + r \\log r)$ algorithm against Garey et al.\\ $O(n \\log n)$.}\n")
        f.write("\\label{tab:benchmark}\n")
        f.write("\\smallskip\n")
        f.write("\\begin{tabular}{llrrrrrr}\n")
        f.write("\\toprule\n")
        f.write("\\textbf{Type} & $n$ & $r$ & $r/n$ & \\textbf{Ours} & \\textbf{Garey} & \\textbf{Speedup} \\\\\n")
        f.write("\\midrule\n")
        
        for ptype in ['convex', 'random_star', 'star']:
            type_results = [r for r in all_results if r['type'] == ptype]
            for r in type_results:
                if r['n'] in [1000, 5000, 10000]:
                    f.write(f"{ptype.replace('_', ' ').title()} & {r['n']:,} & {r['r']:,} & "
                           f"{r['r_ratio']:.0%} & {r['ours_ms']:.2f} & {r['garey_ms']:.2f} & "
                           f"{r['speedup']:.1f}$\\times$ \\\\\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"LaTeX table saved to {latex_path}", file=sys.stderr)


if __name__ == '__main__':
    main()

