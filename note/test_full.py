import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')
from triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)
tris = pt.triangulate()

print('Triangles:', tris)
print('Count:', len(tris))

poly_area = 0
for i in range(len(pts)):
    j = (i + 1) % len(pts)
    poly_area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
poly_area = abs(poly_area) / 2
print('Polygon area:', poly_area)

tri_area = 0
for t in tris:
    a, b, c = pt.pts[t[0]], pt.pts[t[1]], pt.pts[t[2]]
    area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
    tri_area += area
    print(f'  Triangle {t}: area = {area}')
print('Total triangle area:', tri_area)
