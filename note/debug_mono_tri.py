import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import PolygonTriangulator

pts = [(0, 0), (2, 1), (4, 0), (4, 3), (2, 2), (0, 3)]
pt = PolygonTriangulator(pts)

face = [1, 2, 3, 4]
print(f'Face: {face}')
print(f'Face vertices: {[pts[i] for i in face]}')

tris = pt._triangulate_monotone(face)
print(f'Triangles: {tris}')

m = len(face)
indices = face
pts_face = [pts[i] for i in indices]
print(f'Face points: {pts_face}')

sorted_idx = sorted(range(m), key=lambda i: (-pts_face[i][1], pts_face[i][0]))
print(f'Sorted indices (in face): {sorted_idx}')
print(f'Sorted vertices: {[pts_face[i] for i in sorted_idx]}')

top = sorted_idx[0]
bot = sorted_idx[-1]
print(f'Top (in face): {top} = vertex {indices[top]} at {pts_face[top]}')
print(f'Bot (in face): {bot} = vertex {indices[bot]} at {pts_face[bot]}')

left_chain = set()
curr = top
while curr != bot:
    left_chain.add(curr)
    curr = (curr + 1) % m
left_chain.add(bot)
print(f'Left chain (in face): {left_chain}')
