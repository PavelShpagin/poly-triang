from triangulation import *

print('=== Arrow Debug ===')
pts = arrow_shape()
print('Polygon:', pts)

tri = Triangulator(pts)
tri._classify_vertices()
print('Types:', [(i, tri.types[i].name) for i in range(len(pts))])
print('Extrema:', tri.extrema)
print('r:', tri.r)

tri._build_chains()
tri._decompose()
print('Diagonals:', tri.diagonals)

triangles = tri._triangulate_faces()
print('Triangles:', len(triangles), triangles)

print()
print('=== Comb-3 Debug ===')
pts = comb_polygon(3)
print('Polygon:', pts)

tri = Triangulator(pts)
tri._classify_vertices()
print('Types:', [(i, tri.types[i].name) for i in range(len(pts))])
print('Extrema:', tri.extrema)
print('r:', tri.r)

tri._build_chains()
tri._decompose()
print('Diagonals:', tri.diagonals)

triangles = tri._triangulate_faces()
print('Triangles:', len(triangles), triangles)
