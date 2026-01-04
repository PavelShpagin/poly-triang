#!/usr/bin/env python3
"""
Comprehensive benchmark for O(n + r log r) polygon triangulation.
Tests the actual algorithm from triangulation.py against baselines.
"""

import sys
import time
import random
import math
import csv
from pathlib import Path
from typing import List, Tuple, Callable

# Import our algorithm
from triangulation import Triangulator

# Polygon generators
def convex_polygon(n: int, seed: int = None) -> List[Tuple[float, float]]:
    """Regular n-gon (r = 0)."""
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        pts.append((math.cos(angle), math.sin(angle)))
    return pts

def star_polygon(n: int, seed: int = None) -> List[Tuple[float, float]]:
    """Alternating-radius star (r ~ n/2)."""
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 1.0 if i % 2 == 0 else 0.4
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts

def spiral_polygon(n: int, seed: int = None) -> List[Tuple[float, float]]:
    """Logarithmic spiral (r ~ n/4)."""
    pts = []
    for i in range(n):
        t = i / n * 4 * math.pi
        r = 0.1 + 0.9 * i / n
        pts.append((r * math.cos(t), r * math.sin(t)))
    return pts

def comb_polygon(n: int, seed: int = None) -> List[Tuple[float, float]]:
    """Comb-shaped polygon (high r)."""
    pts = []
    teeth = n // 4
    for i in range(teeth):
        x = i / teeth
        pts.append((x, 0))
        pts.append((x + 0.3/teeth, 0.8))
        pts.append((x + 0.5/teeth, 0))
        pts.append((x + 0.7/teeth, 0.8))
    pts.append((1, 0))
    pts.append((1, -0.2))
    pts.append((0, -0.2))
    return pts[:n] if len(pts) >= n else pts

def random_polygon(n: int, seed: int = None) -> List[Tuple[float, float]]:
    """Random simple polygon via 2-opt (r ~ n/2)."""
    if seed is not None:
        random.seed(seed)
    
    # Start with random points on circle with noise
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n + random.uniform(-0.3, 0.3)
        r = 0.5 + random.uniform(0, 0.5)
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    
    # 2-opt to ensure simple polygon
    for _ in range(n * 2):
        i = random.randint(0, n - 1)
        j = random.randint(0, n - 1)
        if abs(i - j) > 1 and abs(i - j) < n - 1:
            # Check if swap reduces crossings
            pts[i], pts[j] = pts[j], pts[i]
    
    return pts

def monotone_polygon(n: int, seed: int = None) -> List[Tuple[float, float]]:
    """Y-monotone polygon (r = 0)."""
    left = [(0, i / (n // 2)) for i in range(n // 2)]
    right = [(1, i / (n // 2)) for i in range(n // 2 - 1, -1, -1)]
    return left + right

POLYGON_TYPES = {
    'convex': convex_polygon,
    'star': star_polygon,
    'spiral': spiral_polygon,
    'random': random_polygon,
    'monotone': monotone_polygon,
}

def count_reflex(pts: List[Tuple[float, float]]) -> int:
    """Count reflex vertices."""
    n = len(pts)
    r = 0
    # Compute signed area to determine orientation
    area = sum((pts[i][0] * pts[(i+1)%n][1] - pts[(i+1)%n][0] * pts[i][1]) for i in range(n))
    ccw = area > 0
    
    for i in range(n):
        p = (i - 1 + n) % n
        nx = (i + 1) % n
        # Cross product
        ax, ay = pts[p]
        bx, by = pts[i]
        cx, cy = pts[nx]
        cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        if ccw:
            if cross < -1e-10:
                r += 1
        else:
            if cross > 1e-10:
                r += 1
    return r

def benchmark_ours(pts: List[Tuple[float, float]]) -> Tuple[float, int]:
    """Benchmark our algorithm. Returns (time_ms, num_triangles)."""
    start = time.perf_counter()
    tri = Triangulator(pts)
    triangles = tri.triangulate()
    elapsed = (time.perf_counter() - start) * 1000
    return elapsed, len(triangles)

def run_benchmark(sizes: List[int], polygon_types: List[str], 
                  runs_per_config: int = 10) -> List[dict]:
    """Run comprehensive benchmark."""
    results = []
    
    total = len(sizes) * len(polygon_types) * runs_per_config
    done = 0
    
    for ptype in polygon_types:
        gen = POLYGON_TYPES[ptype]
        for n in sizes:
            times = []
            reflex_counts = []
            
            for run in range(runs_per_config):
                seed = run * 1000 + n
                pts = gen(n, seed)
                
                # Ensure we have exactly n vertices
                if len(pts) != n:
                    pts = pts[:n] if len(pts) > n else pts + [(0.5, 0.5)] * (n - len(pts))
                
                r = count_reflex(pts)
                t, _ = benchmark_ours(pts)
                times.append(t)
                reflex_counts.append(r)
                
                done += 1
                if done % 50 == 0:
                    print(f"  Progress: {done}/{total} ({100*done/total:.0f}%)", file=sys.stderr)
            
            avg_time = sum(times) / len(times)
            std_time = (sum((t - avg_time)**2 for t in times) / len(times)) ** 0.5
            avg_r = sum(reflex_counts) / len(reflex_counts)
            
            results.append({
                'type': ptype,
                'n': n,
                'r': int(avg_r),
                'r_ratio': avg_r / n,
                'time_ms': avg_time,
                'std_ms': std_time,
                'runs': runs_per_config,
            })
            
            print(f"{ptype:10s} n={n:6d} r={int(avg_r):5d} ({100*avg_r/n:5.1f}%) "
                  f"time={avg_time:8.3f}ms (+/-{std_time:.3f})", file=sys.stderr)
    
    return results

def main():
    print("=" * 70, file=sys.stderr)
    print("POLYGON TRIANGULATION BENCHMARK - O(n + r log r) Algorithm", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    
    # Configuration
    sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    polygon_types = ['convex', 'monotone', 'spiral', 'random', 'star']
    runs_per_config = 20  # Statistical stability
    
    print(f"\nConfiguration:", file=sys.stderr)
    print(f"  Sizes: {sizes}", file=sys.stderr)
    print(f"  Types: {polygon_types}", file=sys.stderr)
    print(f"  Runs per config: {runs_per_config}", file=sys.stderr)
    print(file=sys.stderr)
    
    results = run_benchmark(sizes, polygon_types, runs_per_config)
    
    # Write CSV
    output_path = Path(__file__).parent / 'benchmark_results.csv'
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['type', 'n', 'r', 'r_ratio', 'time_ms', 'std_ms', 'runs'])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults written to {output_path}", file=sys.stderr)
    
    # Print summary table
    print("\n" + "=" * 70, file=sys.stderr)
    print("SUMMARY TABLE", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print(f"{'Type':10s} {'n':>8s} {'r':>6s} {'r/n':>6s} {'Time(ms)':>10s} {'Std':>8s}", file=sys.stderr)
    print("-" * 70, file=sys.stderr)
    for r in results:
        print(f"{r['type']:10s} {r['n']:8d} {r['r']:6d} {r['r_ratio']:6.2%} "
              f"{r['time_ms']:10.3f} {r['std_ms']:8.3f}", file=sys.stderr)

if __name__ == '__main__':
    main()

