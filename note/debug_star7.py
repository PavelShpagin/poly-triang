import sys
import math
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import PolygonTriangulator, create_star_polygon

pts = create_star_polygon(7)
print(f'7-point star: {pts}')
print(f'n = {len(pts)}')

pt = PolygonTriangulator(pts)
tris = pt.triangulate()

# Compute polygon area
n = len(pts)
poly_area = 0
for i in range(n):
    j = (i + 1) % n
    poly_area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
poly_area = abs(poly_area) / 2
print(f'Polygon area: {poly_area}')

# Compute total triangle area
tri_area = 0
for t in tris:
    a, b, c = pt.pts[t[0]], pt.pts[t[1]], pt.pts[t[2]]
    area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
    tri_area += area

print(f'Total triangle area: {tri_area}')
print(f'Triangles: {len(tris)}, expected: {n-2}')
print(f'Area match: {abs(poly_area - tri_area) < 1e-6}')
