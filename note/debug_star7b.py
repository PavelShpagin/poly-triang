import sys
import math
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import PolygonTriangulator, create_star_polygon, VType

pts = create_star_polygon(7)
print(f'n = {len(pts)}')

pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._build_chains()

print(f'Vertex types:')
for i in range(pt.n):
    print(f'  {i}: {pt.V[i].name} at {pts[i]}')

print(f'\nExtrema: {pt.extrema}')
print(f'r = {pt.r}')

print(f'\nChains:')
for cid, chain in pt.chains.items():
    print(f'  Chain {cid}: {chain.vertex_indices}')

pt._make_monotone_optimal()
print(f'\nDiagonals: {pt.diagonals}')
