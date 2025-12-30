#!/usr/bin/env python3
from triangulation import *
from collections import defaultdict
import math

def debug_polygon(name, pts):
    print(f'\n{"="*60}')
    print(f'{name}')
    print(f'{"="*60}')
    print(f'Points: {pts}')
    
    tri = PolygonTriangulator(pts)
    tri._classify_vertices()
    print(f'\nVertex types:')
    for i in range(len(pts)):
        print(f'  {i}: {tri.V[i].name} at {pts[i]}')
    print(f'\nExtrema: {tri.extrema}')
    print(f'Reflex count r = {tri.r}')
    
    if tri.r == 0:
        print('Convex polygon - no decomposition needed')
        return
    
    tri._build_chains()
    print(f'\nChains ({len(tri.chains)} total):')
    for i, c in enumerate(tri.chains):
        print(f'  {i}: {c.vertices}, left_boundary={c.is_left_boundary}, upper={c.upper}, lower={c.lower}')
    
    print(f'\nmin_to_left_chain: {[(k, v.vertices) for k,v in tri.min_to_left_chain.items()]}')
    print(f'max_to_left_chain: {[(k, v.vertices) for k,v in tri.max_to_left_chain.items()]}')
    
    tri._make_monotone()
    print(f'\nDiagonals: {tri.diagonals}')
    
    # Debug face extraction
    adj = defaultdict(set)
    for i in range(tri.n):
        j = tri._nxt(i)
        adj[i].add(j)
        adj[j].add(i)
    for u, v in tri.diagonals:
        adj[u].add(v)
        adj[v].add(u)
    
    print(f'\nAdjacency:')
    for k in sorted(adj.keys()):
        print(f'  {k}: {sorted(adj[k])}')
    
    # Run full triangulation
    triangles = tri.triangulate()
    print(f'\nTriangles ({len(triangles)}):')
    for t in triangles:
        print(f'  {t}')
    
    valid, msg = validate_triangulation(tri.pts, triangles)
    print(f'\nValidation: {msg}')


# Test Arrow
debug_polygon('Arrow', create_arrow_polygon())

# Test 3-tooth comb
debug_polygon('3-tooth Comb', create_comb_polygon(3))

