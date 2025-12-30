from triangulation import PolygonTriangulator, create_paper_example
import math

pts = create_paper_example()
pt = PolygonTriangulator(pts)

# v9 is at (2.5, 4.0)
v9 = pt.pts[9]
print(f'v9 = {v9}')

neighbors = [11, 7, 8, 10]
for n in neighbors:
    vn = pt.pts[n]
    angle = math.atan2(vn[1] - v9[1], vn[0] - v9[0])
    print(f'  v{n} = {vn}, angle = {math.degrees(angle):.1f} deg')

print()
print('Sorted by angle:')
sorted_neighbors = sorted(neighbors, key=lambda x: math.atan2(pt.pts[x][1] - v9[1], pt.pts[x][0] - v9[0]))
print(f'  {sorted_neighbors}')
