#!/usr/bin/env python3
import sys

def ccw(A, B, C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def intersect(A, B, C, D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

# Read polygon
with open(sys.argv[1]) as f:
    n = int(f.readline())
    pts = []
    for _ in range(n):
        x, y = map(float, f.readline().split())
        pts.append((x, y))

# Check for self-intersections
count = 0
for i in range(n):
    for j in range(i+2, n):
        if i == 0 and j == n-1:
            continue
        A, B = pts[i], pts[(i+1)%n]
        C, D = pts[j], pts[(j+1)%n]
        if intersect(A, B, C, D):
            count += 1
            if count < 5:
                print(f'Intersection: edge {i}-{(i+1)%n} with {j}-{(j+1)%n}')

print(f'Total intersections: {count}')
print(f'Polygon vertices: {n}')

