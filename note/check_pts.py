from triangulation import PolygonTriangulator, create_paper_example
pts = create_paper_example()
pt = PolygonTriangulator(pts)
print(pts == pt.pts)
