from triangulation import PolygonTriangulator, create_paper_example
from collections import defaultdict
import math

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._make_monotone()

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
        pt.pts[x][1] - pt.pts[u][1],
        pt.pts[x][0] - pt.pts[u][0]
    ))

# Trace face starting from edge (2, 3)
print('Tracing face starting from edge (2, 3):')
face = [2]
curr = 3
prev = 2
while curr != 2:
    face.append(curr)
    neighbors = adj[curr]
    idx = neighbors.index(prev)
    next_v = neighbors[(idx - 1) % len(neighbors)]
    print(f'  At v{curr}, prev=v{prev}, neighbors={neighbors}, idx={idx}, next=v{next_v}')
    prev = curr
    curr = next_v
    if len(face) > 20:
        print('  Breaking - too long')
        break

print(f'Face: {face}')
