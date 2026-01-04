#!/usr/bin/env python3
"""Quick benchmark with smaller sizes to iterate faster."""
import sys
import time
import math
import random
from triangulate_correct import triangulate, gen_convex, gen_star, gen_random_star, gen_spiral

def bench(gen, n, runs=10):
    times = []
    r_vals = []
    for seed in range(runs):
        pts = gen(n, seed)
        t0 = time.perf_counter()
        _, r = triangulate(pts)
        times.append((time.perf_counter() - t0) * 1000)
        r_vals.append(r)
    return sum(times)/len(times), sum(r_vals)/len(r_vals)

print("Quick benchmark (10 runs each):")
print(f"{'Type':<12} {'n':>6} {'r':>6} {'time_ms':>10}")
print("-" * 40)

for ptype, gen in [('convex', gen_convex), ('spiral', gen_spiral), 
                   ('random_star', gen_random_star), ('star', gen_star)]:
    for n in [100, 500, 1000, 2000, 5000, 10000]:
        t, r = bench(gen, n, runs=10)
        print(f"{ptype:<12} {n:>6} {int(r):>6} {t:>10.3f}")
    print()

