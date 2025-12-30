from triangulation import PolygonTriangulator, create_paper_example

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()
pt._make_monotone()

# The problematic face (using pt.pts indices)
face = [2, 3, 4, 5, 6, 7, 9, 11]
indices = face
m = len(indices)

print('Face:', face)
print('Face vertices:')
for i, idx in enumerate(face):
    print(f'  face[{i}] = v{idx} = {pt.pts[idx]}')

local_pts = [pt.pts[i] for i in indices]
sorted_idx = sorted(range(m), key=lambda i: (-local_pts[i][1], local_pts[i][0]))

print()
print('Sorted order (by y descending):')
for i, si in enumerate(sorted_idx):
    print(f'  {i}: face[{si}] = v{indices[si]} = {local_pts[si]}')

top = sorted_idx[0]
bot = sorted_idx[-1]
print(f'Top: face[{top}] = v{indices[top]}')
print(f'Bot: face[{bot}] = v{indices[bot]}')

# Build chains
left_chain = set()
curr = top
while curr != bot:
    left_chain.add(curr)
    curr = (curr + 1) % m
left_chain.add(bot)

right_chain = set()
curr = top
while curr != bot:
    right_chain.add(curr)
    curr = (curr - 1 + m) % m
right_chain.add(bot)

print()
print('Left chain (face indices):', sorted(left_chain))
print('Right chain (face indices):', sorted(right_chain))
