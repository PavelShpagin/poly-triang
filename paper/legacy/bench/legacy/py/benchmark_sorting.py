#!/usr/bin/env python3
"""
Benchmark specifically comparing sorting costs:
- Our algorithm: Sort only r reflex vertices
- Garey: Sort all n vertices

This isolates the O(r log r) vs O(n log n) advantage.
"""
import sys
import time
import random
import math
from triangulate_correct import gen_convex, gen_star, gen_random_star
from triangulate_correct import classify_vertices, signed_area


def time_garey_sort(pts):
    """Time sorting ALL vertices (Garey approach)."""
    n = len(pts)
    t0 = time.perf_counter()
    sorted_verts = sorted(range(n), key=lambda i: (-pts[i][1], pts[i][0]))
    return (time.perf_counter() - t0) * 1000, n


def time_our_sort(pts):
    """Time sorting ONLY reflex vertices (our approach)."""
    from triangulate_correct import VertexType
    
    if signed_area(pts) < 0:
        pts = list(reversed(pts))
    
    n = len(pts)
    t0 = time.perf_counter()
    types, r = classify_vertices(pts)
    
    # Get extrema (reflex + start/end vertices)
    extrema = [i for i in range(n) if types[i] in 
               (VertexType.START, VertexType.END, VertexType.SPLIT, VertexType.MERGE)]
    
    # Sort only extrema
    sorted_extrema = sorted(extrema, key=lambda i: (-pts[i][1], pts[i][0]))
    elapsed = (time.perf_counter() - t0) * 1000
    
    return elapsed, len(extrema)


def run_benchmark():
    """Run sorting benchmark."""
    print("=" * 80)
    print("SORTING COST COMPARISON: O(r log r) vs O(n log n)")
    print("=" * 80)
    
    sizes = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
    
    for ptype, gen in [('convex', gen_convex), ('random_star', gen_random_star), ('star', gen_star)]:
        print(f"\n{ptype.upper()}:")
        print(f"{'n':>8} {'r':>8} {'Ours(ms)':>12} {'Garey(ms)':>12} {'Speedup':>10} {'Theory':>10}")
        print("-" * 70)
        
        for n in sizes:
            pts = gen(n, seed=42)
            
            # Run multiple times for stability
            ours_times = []
            garey_times = []
            r_val = 0
            
            for _ in range(5):
                t_ours, r_val = time_our_sort(pts)
                t_garey, _ = time_garey_sort(pts)
                ours_times.append(t_ours)
                garey_times.append(t_garey)
            
            avg_ours = sum(ours_times) / len(ours_times)
            avg_garey = sum(garey_times) / len(garey_times)
            
            speedup = avg_garey / avg_ours if avg_ours > 0 else float('inf')
            
            # Theoretical speedup: (n log n) / (r log r) for sorting
            if r_val > 1:
                theory = (n * math.log2(n)) / (r_val * math.log2(r_val))
            else:
                theory = n  # r=0 means we skip entirely
            
            print(f"{n:>8} {r_val:>8} {avg_ours:>12.3f} {avg_garey:>12.3f} "
                  f"{speedup:>10.2f}x {theory:>10.1f}x")
    
    print("\n" + "=" * 80)
    print("Note: 'Ours' includes classification + sorting extrema")
    print("      'Garey' is sorting all n vertices")
    print("      For r << n, our sorting cost is significantly lower")
    print("=" * 80)


if __name__ == '__main__':
    run_benchmark()

