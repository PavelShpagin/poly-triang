"""Debug face extraction."""

import sys
import math
from collections import defaultdict

sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import PolygonTriangulator

def create_monotone_polygon():
    return [
        (0, 0), (2, 1), (4, 0), (4, 3), (2, 2), (0, 3)
    ]

pts = create_monotone_polygon()
n = len(pts)

pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._build_chains()
pt._make_monotone_optimal()

print(f"Diagonals: {pt.diagonals}")

# Build adjacency
adj = defaultdict(list)
for i in range(n):
    j = (i + 1) % n
    adj[i].append(j)
    adj[j].append(i)
    print(f"Edge {i}-{j}")

for u, v in pt.diagonals:
    adj[u].append(v)
    adj[v].append(u)
    print(f"Diagonal {u}-{v}")

print("\nAdjacency before sorting:")
for u in adj:
    print(f"  {u}: {adj[u]}")

# Sort neighbors by angle
for u in adj:
    adj[u].sort(key=lambda x: math.atan2(
        pts[x][1] - pts[u][1],
        pts[x][0] - pts[u][0]
    ))

print("\nAdjacency after sorting:")
for u in adj:
    angles = [(x, math.atan2(pts[x][1] - pts[u][1], pts[x][0] - pts[u][0])) for x in adj[u]]
    print(f"  {u}: {angles}")

# Extract faces
used = set()
faces = []

for i in range(n):
    for start_v in adj[i]:
        if (i, start_v) in used:
            continue
        print(f"\nStarting face from edge ({i}, {start_v})")
        face = [i]
        curr = start_v
        prev = i
        steps = 0
        while curr != i and steps < 20:
            face.append(curr)
            used.add((prev, curr))
            neighbors = adj[curr]
            idx = neighbors.index(prev)
            next_v = neighbors[(idx - 1) % len(neighbors)]
            print(f"  At {curr}, prev={prev}, neighbors={neighbors}, idx={idx}, next={next_v}")
            prev = curr
            curr = next_v
            steps += 1
        used.add((prev, curr))
        print(f"  Face: {face}")
        
        if len(face) >= 3:
            pts_face = [pts[k] for k in face]
            area = 0
            for j in range(len(pts_face)):
                k = (j + 1) % len(pts_face)
                area += pts_face[j][0] * pts_face[k][1] - pts_face[k][0] * pts_face[j][1]
            print(f"  Area: {area}")
            if area > 1e-9:
                faces.append(face)

print(f"\nFinal faces: {faces}")

