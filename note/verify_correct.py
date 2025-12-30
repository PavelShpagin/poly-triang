from triangulation import PolygonTriangulator, create_paper_example, validate_triangulation

pts = create_paper_example()
pt = PolygonTriangulator(pts)
tris = pt.triangulate()

print('Using pt.pts (CCW version):')
valid, msg = validate_triangulation(pt.pts, tris)
print(f'Valid: {valid}, Msg: {msg}')

# Compute areas using pt.pts
poly_area = 0
for i in range(len(pt.pts)):
    j = (i + 1) % len(pt.pts)
    poly_area += pt.pts[i][0] * pt.pts[j][1] - pt.pts[j][0] * pt.pts[i][1]
poly_area = abs(poly_area) / 2
print(f'Polygon area: {poly_area}')

tri_area = 0
for t in tris:
    a, b, c = pt.pts[t[0]], pt.pts[t[1]], pt.pts[t[2]]
    area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
    tri_area += area
print(f'Triangle area: {tri_area}')
