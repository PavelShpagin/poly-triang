import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')
from triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)

face = [2, 11, 9, 7, 6, 5, 4, 3]
print('Face:', face)
print('Face vertices:', [pts[i] for i in face])

tris = pt._triangulate_monotone(face)
print('Triangles:', tris)

face_pts = [pts[i] for i in face]
face_area = 0
for j in range(len(face_pts)):
    k = (j + 1) % len(face_pts)
    face_area += face_pts[j][0] * face_pts[k][1] - face_pts[k][0] * face_pts[j][1]
face_area = abs(face_area) / 2
print('Face area:', face_area)

tri_area = 0
for t in tris:
    a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
    area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
    tri_area += area
    print('  Triangle', t, ': area =', area)
print('Total triangle area:', tri_area)
