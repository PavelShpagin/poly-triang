#!/usr/bin/env python3
"""Quick test of the triangulation algorithm."""
import math
from triangulation import Triangulator

# Test on a simple star polygon (non-convex)
n = 10
pts = []
for i in range(n):
    angle = 2 * math.pi * i / n
    r = 1.0 if i % 2 == 0 else 0.4
    pts.append((r * math.cos(angle), r * math.sin(angle)))

tri = Triangulator(pts)
triangles = tri.triangulate()
print(f'Star polygon: n={n}, r={tri.r}, triangles={len(triangles)}')
print(f'Expected triangles: {n-2}')

# Validate: should have n-2 triangles
assert len(triangles) == n - 2, f'Wrong triangle count: {len(triangles)} vs {n-2}'
print('PASS: Star polygon')

# Test random polygon
import random
random.seed(42)
n = 100
pts = []
for i in range(n):
    angle = 2 * math.pi * i / n + random.uniform(-0.2, 0.2)
    radius = 0.5 + random.uniform(0, 0.5)
    pts.append((radius * math.cos(angle), radius * math.sin(angle)))

tri = Triangulator(pts)
triangles = tri.triangulate()
print(f'Random polygon: n={n}, r={tri.r}, triangles={len(triangles)}')
assert len(triangles) == n - 2, f'Wrong triangle count: {len(triangles)} vs {n-2}'
print('PASS: Random polygon')

# Test convex polygon
n = 50
pts = [(math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)) for i in range(n)]
tri = Triangulator(pts)
triangles = tri.triangulate()
print(f'Convex polygon: n={n}, r={tri.r}, triangles={len(triangles)}')
assert len(triangles) == n - 2
assert tri.r == 0
print('PASS: Convex polygon')

print('\nAll tests passed!')

