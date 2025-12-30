import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')
from triangulation import PolygonTriangulator, create_paper_example
from collections import defaultdict
import math

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._make_monotone()

print('Diagonals:', pt.diagonals)

# Build adjacency
adj = defaultdict(list)
for i in range(pt.n):
    j = (i + 1) % pt.n
    adj[i].append(j)
    adj[j].append(i)
for u, v in pt.diagonals:
    adj[u].append(v)
    adj[v].append(u)

for u in adj:
    adj[u].sort(key=lambda x: math.atan2(
        pts[x][1] - pts[u][1],
        pts[x][0] - pts[u][0]
    ))

print('Adjacency:')
for u in sorted(adj.keys()):
    print(f'  {u}: {adj[u]}')

# Extract faces
used = set()
faces = []

for i in range(pt.n):
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
            if len(face) > 2 * pt.n:
                break
        used.add((prev, curr))
        
        pts_face = [pts[k] for k in face]
        area = 0
        for j in range(len(pts_face)):
            k = (j + 1) % len(pts_face)
            area += pts_face[j][0] * pts_face[k][1] - pts_face[k][0] * pts_face[j][1]
        
        print(f'Face {face}: area = {area/2}')
        if len(face) >= 3 and area > 1e-9:
            faces.append(face)

print('\nFaces to triangulate:', faces)
