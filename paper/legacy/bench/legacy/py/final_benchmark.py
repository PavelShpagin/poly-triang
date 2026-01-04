#!/usr/bin/env python3
"""
Final comprehensive benchmark for the paper.
Tests multiple polygon types with statistical stability.
"""
import subprocess
import sys
import time
import csv
import random
import math
import tempfile
import os
from pathlib import Path

ROOT = Path(__file__).parent.parent
BIN_DIR = ROOT / 'build' / 'bin'
RESULTS_DIR = ROOT / 'results'

def gen_convex(n, seed=None):
    """Generate convex polygon with n vertices."""
    if seed is not None:
        random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100 + random.uniform(-1, 1)  # Slight perturbation
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts

def gen_star(n, seed=None):
    """Generate star-shaped polygon (alternating radii)."""
    if seed is not None:
        random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100 if i % 2 == 0 else 40
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts

def gen_spiral(n, seed=None):
    """Generate monotone spiral polygon."""
    if seed is not None:
        random.seed(seed)
    pts = []
    for i in range(n):
        t = i / n
        angle = 4 * math.pi * t
        r = 20 + 80 * t
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts

def gen_random_star(n, seed=None):
    """Generate star polygon with random radii."""
    if seed is not None:
        random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = random.uniform(30, 100)
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts

def gen_comb(n, seed=None):
    """Generate comb-shaped polygon."""
    if seed is not None:
        random.seed(seed)
    pts = []
    teeth = n // 4
    for i in range(teeth):
        x = i * 10
        pts.append((x, 0))
        pts.append((x + 2, 50))
        pts.append((x + 5, 50))
        pts.append((x + 7, 0))
    return pts[:n]

def write_polygon(pts, path):
    """Write polygon to file."""
    with open(path, 'w') as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")

def run_algo(algo, input_path, output_path):
    """Run algorithm and return time in ms."""
    cli = BIN_DIR / f"{algo}_cli"
    if not cli.exists():
        return None, None
    
    try:
        result = subprocess.run(
            [str(cli), '--input', str(input_path), '--output', str(output_path)],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0:
            output = result.stdout
            # Parse time
            for line in output.split('\n'):
                if 'time_ms=' in line:
                    time_ms = float(line.split('time_ms=')[1].split(',')[0].split()[0])
                    # Parse reflex count if available
                    r = 0
                    if 'reflex_count=' in line:
                        r = int(line.split('reflex_count=')[1].split(',')[0].split()[0])
                    return time_ms, r
    except subprocess.TimeoutExpired:
        pass
    return None, None

def main():
    print("=" * 80)
    print("FINAL COMPREHENSIVE BENCHMARK")
    print("=" * 80)
    
    RESULTS_DIR.mkdir(exist_ok=True)
    
    # Configuration
    generators = {
        'convex': gen_convex,
        'spiral': gen_spiral,
        'random_star': gen_random_star,
        'star': gen_star,
        'comb': gen_comb,
    }
    
    sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    polygons_per_config = 50  # Statistical stability
    algorithms = ['reflex', 'earcut']
    
    # Results storage
    csv_path = RESULTS_DIR / 'final_benchmark.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['polygon_type', 'n', 'algorithm', 'run', 'time_ms', 'reflex_count'])
    
    # Run benchmarks
    with tempfile.TemporaryDirectory() as tmpdir:
        for ptype, gen in generators.items():
            print(f"\n{ptype.upper()}:")
            
            for n in sizes:
                reflex_times = []
                earcut_times = []
                r_vals = []
                
                for run in range(polygons_per_config):
                    pts = gen(n, seed=run)
                    input_path = Path(tmpdir) / f"poly_{ptype}_{n}_{run}.poly"
                    output_path = Path(tmpdir) / f"out_{ptype}_{n}_{run}.tri"
                    
                    write_polygon(pts, input_path)
                    
                    for algo in algorithms:
                        time_ms, r = run_algo(algo, input_path, output_path)
                        
                        if time_ms is not None:
                            with open(csv_path, 'a', newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow([ptype, n, algo, run, time_ms, r or 0])
                            
                            if algo == 'reflex':
                                reflex_times.append(time_ms)
                                if r is not None:
                                    r_vals.append(r)
                            elif algo == 'earcut':
                                earcut_times.append(time_ms)
                
                if reflex_times and earcut_times:
                    avg_reflex = sum(reflex_times) / len(reflex_times)
                    avg_earcut = sum(earcut_times) / len(earcut_times)
                    avg_r = sum(r_vals) / len(r_vals) if r_vals else 0
                    speedup = avg_earcut / avg_reflex if avg_reflex > 0 else 0
                    
                    print(f"  n={n:>6}  r={int(avg_r):>5}  "
                          f"Ours={avg_reflex:>8.3f}ms  Earcut={avg_earcut:>8.3f}ms  "
                          f"Speedup={speedup:>6.2f}x  ({len(reflex_times)} samples)")
    
    print(f"\nResults saved to {csv_path}")
    
    # Generate summary
    import pandas as pd
    df = pd.read_csv(csv_path)
    
    summary = df.groupby(['polygon_type', 'n', 'algorithm']).agg({
        'time_ms': ['mean', 'std'],
        'reflex_count': 'first'
    }).reset_index()
    summary.columns = ['polygon_type', 'n', 'algorithm', 'time_mean', 'time_std', 'reflex_count']
    
    summary_path = RESULTS_DIR / 'final_summary.csv'
    summary.to_csv(summary_path, index=False)
    print(f"Summary saved to {summary_path}")

if __name__ == '__main__':
    main()

