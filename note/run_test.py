#!/usr/bin/env python3
import sys
import math
from collections import defaultdict
sys.path.insert(0, '.')
from note.triangulation import PolygonTriangulator, create_paper_example, validate_triangulation, VType

pts = create_paper_example()
pt = PolygonTriangulator(pts)

# Show classification
pt._classify_vertices()
print("=== Vertex Classification ===")
for v in pt.V:
    print(f"  v{v.idx}: ({v.x:.1f}, {v.y:.1f}) {v.vtype.name}")

print()
print("Reflex vertices:")
for v in pt.V:
    if v.vtype in (VType.SPLIT, VType.MERGE):
        print(f"  v{v.idx}: {v.vtype.name}")

# Run monotone decomposition
pt2 = PolygonTriangulator(pts)
pt2._classify_vertices()
pt2._make_monotone()

print()
print(f"=== Diagonals ===")
print(f"Diagonals: {pt2.diagonals}")

# Debug face extraction
n = pt2.n
adj = defaultdict(list)
for i in range(n):
    j = (i + 1) % n
    adj[i].append(j)
    adj[j].append(i)
for u, v in pt2.diagonals:
    adj[u].append(v)
    adj[v].append(u)

for u in adj:
    adj[u].sort(key=lambda x: math.atan2(
        pt2.pts[x][1] - pt2.pts[u][1],
        pt2.pts[x][0] - pt2.pts[u][0]
    ))

used = set()
faces = []

print()
print("=== Face Extraction ===")
for i in range(n):
    for start_v in adj[i]:
        if (i, start_v) in used:
            continue
        face = [i]
        curr = start_v
        prev = i
        while curr != i:
            face.append(curr)
            used.add((prev, curr))
            neighbors = adj[curr]
            idx = neighbors.index(prev)
            next_v = neighbors[(idx - 1) % len(neighbors)]
            prev = curr
            curr = next_v
            if len(face) > 2 * n:
                break
        used.add((prev, curr))
        
        if len(face) >= 3:
            pts_face = [pt2.pts[k] for k in face]
            area = 0
            for j in range(len(pts_face)):
                k = (j + 1) % len(pts_face)
                area += pts_face[j][0] * pts_face[k][1] - pts_face[k][0] * pts_face[j][1]
            if area > 1e-9:
                faces.append(face)
                print(f"  Interior: {face} (area={area:.2f})")
            else:
                print(f"  Exterior: {face} (area={area:.2f})")

print(f"\n{len(faces)} interior faces")

# Run full triangulation
pt3 = PolygonTriangulator(pts)
tris = pt3.triangulate()
valid, msg = validate_triangulation(pt3.pts, tris)

print()
print(f"=== Results ===")
print(f"Valid: {valid}, Msg: {msg}")
print(f"Triangles: {len(tris)} (expected {len(pts) - 2})")

