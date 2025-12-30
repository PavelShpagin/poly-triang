#!/usr/bin/env python3
"""
Quick benchmark to verify O(n + r log r) behavior.
Profiles individual phases of the algorithm.
"""

import sys
import time
import math
import random
from typing import List, Tuple

from triangulate_correct import (
    classify_vertices, make_monotone, extract_faces, 
    triangulate_monotone, signed_area, gen_convex, gen_star, gen_random_star
)


def profile_phases(pts: List[Tuple[float, float]]) -> dict:
    """Profile each phase of the algorithm."""
    n = len(pts)
    if signed_area(pts) < 0:
        pts = list(reversed(pts))
    
    # Phase 1: Classify vertices - O(n)
    t0 = time.perf_counter()
    types, r = classify_vertices(pts)
    t_classify = (time.perf_counter() - t0) * 1000
    
    if r == 0:
        return {'n': n, 'r': 0, 'classify_ms': t_classify, 'decompose_ms': 0, 
                'extract_ms': 0, 'triangulate_ms': 0, 'total_ms': t_classify}
    
    # Phase 2: Monotone decomposition - O(n + r log r)
    t0 = time.perf_counter()
    diagonals = make_monotone(pts, types)
    t_decompose = (time.perf_counter() - t0) * 1000
    
    # Phase 3: Extract faces - should be O(n + diagonals)
    t0 = time.perf_counter()
    faces = extract_faces(pts, diagonals)
    t_extract = (time.perf_counter() - t0) * 1000
    
    # Phase 4: Triangulate faces - O(n)
    t0 = time.perf_counter()
    triangles = []
    for face in faces:
        triangles.extend(triangulate_monotone(pts, face))
    t_tri = (time.perf_counter() - t0) * 1000
    
    return {
        'n': n, 'r': r, 'diagonals': len(diagonals), 'faces': len(faces),
        'triangles': len(triangles),
        'classify_ms': t_classify,
        'decompose_ms': t_decompose,
        'extract_ms': t_extract,
        'triangulate_ms': t_tri,
        'total_ms': t_classify + t_decompose + t_extract + t_tri
    }


def main():
    print("=" * 80)
    print("PHASE PROFILING - O(n + r log r) Algorithm")
    print("=" * 80)
    
    sizes = [100, 500, 1000, 2000, 5000, 10000]
    
    for ptype, gen in [('convex', gen_convex), ('star', gen_star), ('random_star', gen_random_star)]:
        print(f"\n{ptype.upper()}:")
        print(f"{'n':>6} {'r':>6} {'classify':>10} {'decompose':>10} {'extract':>10} {'triangulate':>12} {'total':>10}")
        print("-" * 70)
        
        for n in sizes:
            pts = gen(n, seed=42)
            result = profile_phases(pts)
            
            print(f"{result['n']:>6} {result['r']:>6} "
                  f"{result['classify_ms']:>10.3f} {result['decompose_ms']:>10.3f} "
                  f"{result['extract_ms']:>10.3f} {result['triangulate_ms']:>12.3f} "
                  f"{result['total_ms']:>10.3f}")
    
    print("\n" + "=" * 80)
    print("SCALING ANALYSIS")
    print("=" * 80)
    
    # Compare scaling between convex (r=0) and star (r~n/2)
    print("\nTime ratio (star/convex) - should be O(log n) if algorithm is correct:")
    for n in [1000, 2000, 5000, 10000]:
        pts_convex = gen_convex(n)
        pts_star = gen_star(n)
        
        t_convex = profile_phases(pts_convex)['total_ms']
        t_star = profile_phases(pts_star)['total_ms']
        
        ratio = t_star / t_convex if t_convex > 0 else float('inf')
        # Expected ratio: (n + r log r) / n ≈ 1 + 0.5 * log(n/2) ≈ 1 + 0.5 * (log n - 1)
        expected = 1 + 0.5 * (math.log2(n) - 1)
        
        print(f"  n={n:>5}: ratio={ratio:.2f}x (expected ~{expected:.2f}x for O(n + r log r))")


if __name__ == '__main__':
    main()

