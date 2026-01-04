#!/usr/bin/env python3
"""
Comprehensive benchmark for O(n + r log r) polygon triangulation.
Uses correct implementation and proper simple polygon generators.
"""

import sys
import time
import csv
import math
import random
from pathlib import Path
from typing import List, Tuple

from triangulate_correct import triangulate, gen_convex, gen_star, gen_random_star, gen_spiral


def benchmark(pts: List[Tuple[float, float]], warmup: int = 1, runs: int = 5) -> Tuple[float, float, int]:
    """Benchmark triangulation. Returns (mean_ms, std_ms, reflex_count)."""
    # Warmup
    for _ in range(warmup):
        triangulate(pts)
    
    times = []
    r = 0
    for _ in range(runs):
        start = time.perf_counter()
        _, r = triangulate(pts)
        elapsed = (time.perf_counter() - start) * 1000
        times.append(elapsed)
    
    mean = sum(times) / len(times)
    std = (sum((t - mean)**2 for t in times) / len(times)) ** 0.5
    return mean, std, r


def main():
    print("=" * 75, file=sys.stderr)
    print("POLYGON TRIANGULATION BENCHMARK - O(n + r log r)", file=sys.stderr)
    print("=" * 75, file=sys.stderr)
    
    # Test configurations
    sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    polygon_types = {
        'convex': gen_convex,
        'spiral': gen_spiral,
        'random_star': gen_random_star,
        'star': gen_star,
    }
    
    runs_per_size = 20  # Multiple polygons per size for statistics
    timing_runs = 5     # Timing runs per polygon
    
    results = []
    
    for ptype, gen in polygon_types.items():
        print(f"\n{ptype.upper()}:", file=sys.stderr)
        print(f"{'n':>8s} {'r':>8s} {'r/n':>8s} {'time_ms':>12s} {'std':>10s}", file=sys.stderr)
        print("-" * 50, file=sys.stderr)
        
        for n in sizes:
            all_times = []
            all_r = []
            
            for seed in range(runs_per_size):
                pts = gen(n, seed=seed)
                mean_ms, std_ms, r = benchmark(pts, warmup=1, runs=timing_runs)
                all_times.append(mean_ms)
                all_r.append(r)
            
            avg_time = sum(all_times) / len(all_times)
            std_time = (sum((t - avg_time)**2 for t in all_times) / len(all_times)) ** 0.5
            avg_r = sum(all_r) / len(all_r)
            
            results.append({
                'type': ptype,
                'n': n,
                'r': int(avg_r),
                'r_ratio': avg_r / n,
                'time_ms': avg_time,
                'std_ms': std_time,
                'runs': runs_per_size * timing_runs,
            })
            
            print(f"{n:>8d} {int(avg_r):>8d} {avg_r/n:>8.2%} {avg_time:>12.4f} {std_time:>10.4f}", file=sys.stderr)
    
    # Write CSV
    output_path = Path(__file__).parent / 'benchmark_final.csv'
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['type', 'n', 'r', 'r_ratio', 'time_ms', 'std_ms', 'runs'])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults saved to {output_path}", file=sys.stderr)
    
    # Print summary
    print("\n" + "=" * 75, file=sys.stderr)
    print("SUMMARY: Time (ms) at n=10000", file=sys.stderr)
    print("=" * 75, file=sys.stderr)
    for r in results:
        if r['n'] == 10000:
            print(f"  {r['type']:15s} r={r['r']:5d} ({r['r_ratio']:5.1%})  time={r['time_ms']:8.3f}ms", file=sys.stderr)


if __name__ == '__main__':
    main()

