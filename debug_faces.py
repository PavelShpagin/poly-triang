import sys
import math
from collections import defaultdict
sys.path.insert(0, .)
from note.triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._make_monotone()

print(Diagonals:, pt.diagonals)

n = pt.n
adj = defaultdict(list)
for i in range(n):
    j = (i + 1) % n
    adj[i].append(j)
    adj[j].append(i)
for u, v in pt.diagonals:
    adj[u].append(v)
    adj[v].append(u)

for u in adj:
    adj[u].sort(key=lambda x: math.atan2(
        pt.pts[x][1] - pt.pts[u][1],
        pt.pts[x][0] - pt.pts[u][0]
    ))

used = set()
faces = []

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
            pts_face = [pt.pts[k] for k in face]
            area = 0
            for j in range(len(pts_face)):
                k = (j + 1) % len(pts_face)
                area += pts_face[j][0] * pts_face[k][1] - pts_face[k][0] * pts_face[j][1]
            if area > 1e-9:
                faces.append(face)
                print(Interior:, face)
            else:
                print(Exterior:, face)

print(len(faces), interior
