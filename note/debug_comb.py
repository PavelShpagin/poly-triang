import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import PolygonTriangulator, create_comb_polygon

pts = create_comb_polygon(3)
print(f'Comb polygon: {pts}')
print(f'n = {len(pts)}')

pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._build_chains()
pt._make_monotone_optimal()

print(f'Extrema: {pt.extrema}')
print(f'r = {pt.r}')
print(f'Diagonals: {pt.diagonals}')

tris = pt._triangulate_monotone_pieces()
print(f'Triangles ({len(tris)}): {tris}')

# Check for degenerate triangles
for t in tris:
    a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
    area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    if abs(area) < 1e-9:
        print(f'DEGENERATE: {t} -> {a}, {b}, {c}')
