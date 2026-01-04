#!/usr/bin/env python3
import subprocess
import time
import random
import math
import sys

n = 10000
seed = 42
poly_file = '/tmp/test_poly.poly'
out_file = '/tmp/test_out.tri'

# Generate test polygon
random.seed(seed)
angles = sorted([random.uniform(0, 2*math.pi) for _ in range(n)])
with open(poly_file, 'w') as f:
    f.write(f'{n}\n')
    for a in angles:
        r = 100 + random.uniform(0, 50)
        x = 500 + r * math.cos(a)
        y = 500 + r * math.sin(a)
        f.write(f'{x} {y}\n')

print(f"Testing with n={n}, seed={seed}")

# Time ours
times_ours = []
for _ in range(10):
    t0 = time.time()
    r = subprocess.run(['./build/bin/reflex_cli', '--input', poly_file, '--output', out_file],
                       capture_output=True, text=True)
    t1 = time.time()
    if r.returncode == 0:
        times_ours.append((t1-t0)*1000)
        
# Time polypartition
times_pp = []
for _ in range(10):
    t0 = time.time()
    r = subprocess.run(['./build/bin/polypartition_mono_cli', '--input', poly_file, '--output', out_file],
                       capture_output=True, text=True)
    t1 = time.time()
    if r.returncode == 0:
        times_pp.append((t1-t0)*1000)

print(f"Ours:          {sum(times_ours)/len(times_ours):.2f} ms avg (10 runs)")
print(f"PolyPartition: {sum(times_pp)/len(times_pp):.2f} ms avg (10 runs)")
print(f"Speedup: {100*(1 - sum(times_ours)/len(times_ours) / (sum(times_pp)/len(times_pp))):.1f}%")

