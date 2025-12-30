from triangulation import PolygonTriangulator, create_paper_example
pts = create_paper_example()
pt = PolygonTriangulator(pts)

def signed_area_2(pts):
    s = 0
    for i in range(len(pts)):
        j = (i + 1) % len(pts)
        s += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
    return s

print('Original area:', signed_area_2(pts))
print('pt.pts area:', signed_area_2(pt.pts))
print('Reversed?', pts[::-1] == pt.pts)
