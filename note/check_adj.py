from triangulation import PolygonTriangulator, create_paper_example
from collections import defaultdict
import math

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._make_monotone()

print('Diagonals:', pt.diagonals)

adj = defaultdict(list)
for i in range(pt.n):
    j = (i + 1) % pt.n
    adj[i].append(j)
    adj[j].append(i)
for u, v in pt.diagonals:
    adj[u].append(v)
    adj[v].append(u)

print()
print('Adjacency before sorting:')
for u in sorted(adj.keys()):
    print(f'  v{u}: {adj[u]}')

for u in adj:
    adj[u].sort(key=lambda x: math.atan2(
        pt.pts[x][1] - pt.pts[u][1],
        pt.pts[x][0] - pt.pts[u][0]
    ))

print()
print('Adjacency after sorting:')
for u in sorted(adj.keys()):
    print(f'  v{u}: {adj[u]}')
