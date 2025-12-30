import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')
from triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)

faces = [[0, 11, 2, 1], [2, 11, 9, 7, 6, 5, 4, 3], [7, 9, 8], [9, 11, 10]]

total_face_area = 0
total_tri_area = 0

for face in faces:
    face_pts = [pts[i] for i in face]
    face_area = 0
    for j in range(len(face_pts)):
        k = (j + 1) % len(face_pts)
        face_area += face_pts[j][0] * face_pts[k][1] - face_pts[k][0] * face_pts[j][1]
    face_area = abs(face_area) / 2
    total_face_area += face_area
    
    tris = pt._triangulate_monotone(face)
    
    tri_area = 0
    for t in tris:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
        tri_area += area
    total_tri_area += tri_area
    
    print(f'Face {face}: face_area={face_area}, tri_area={tri_area}, tris={len(tris)}')

print(f'Total face area: {total_face_area}')
print(f'Total tri area: {total_tri_area}')

# Compare with polygon area
poly_area = 0
for i in range(len(pts)):
    j = (i + 1) % len(pts)
    poly_area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
poly_area = abs(poly_area) / 2
print(f'Polygon area: {poly_area}')
