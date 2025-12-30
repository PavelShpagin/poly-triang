from triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)
tris = pt.triangulate()

print('Triangles:')
total = 0
for t in tris:
    a, b, c = pt.pts[t[0]], pt.pts[t[1]], pt.pts[t[2]]
    area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    signed_area = area / 2
    abs_area = abs(area) / 2
    total += abs_area
    print(f'  {t}: {a} -> {b} -> {c}, signed={signed_area:.2f}, abs={abs_area:.2f}')

print(f'Total: {total}')
