from triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._make_monotone()

print('Diagonals:', pt.diagonals)
print('pt.pts (CCW):')
for i, p in enumerate(pt.pts):
    print(f'  v{i}: {p}, type={pt.V[i].name}')
