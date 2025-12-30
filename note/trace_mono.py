import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

pts = [(0, 0), (2, 1), (4, 0), (4, 3), (2, 2), (0, 3)]

face = [1, 2, 3, 4]
indices = face
m = len(face)
pts_face = [pts[i] for i in indices]

print(f'Face: {face}')
print(f'Face points: {pts_face}')

sorted_idx = sorted(range(m), key=lambda i: (-pts_face[i][1], pts_face[i][0]))
print(f'Sorted indices: {sorted_idx}')

top = sorted_idx[0]
bot = sorted_idx[-1]
print(f'top={top}, bot={bot}')

left_chain = set()
curr = top
while curr != bot:
    left_chain.add(curr)
    curr = (curr + 1) % m
left_chain.add(bot)
print(f'left_chain: {left_chain}')

stack = [sorted_idx[0], sorted_idx[1]]
tris = []
print(f'Initial stack: {stack}')

def triangle_area_2(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

for i in range(2, m):
    u = sorted_idx[i]
    u_left = u in left_chain
    top_left = stack[-1] in left_chain
    print(f'\ni={i}, u={u} (poly vertex {indices[u]}), u_left={u_left}, top_left={top_left}')
    print(f'Stack: {stack}')
    
    if u_left != top_left:
        print('Different chains')
        while len(stack) > 1:
            v = stack.pop()
            tri = (indices[u], indices[v], indices[stack[-1]])
            print(f'  Triangle: {tri}')
            tris.append(tri)
        stack.pop()
        stack.append(sorted_idx[i-1])
        stack.append(u)
    else:
        print('Same chain')
        last = stack.pop()
        while stack:
            v = stack[-1]
            area = triangle_area_2(pts_face[u], pts_face[last], pts_face[v])
            valid = (area > 1e-9) if u_left else (area < -1e-9)
            print(f'  Check: u={u}, last={last}, v={v}, area={area}, valid={valid}')
            if valid:
                tri = (indices[u], indices[last], indices[v])
                print(f'  Triangle: {tri}')
                tris.append(tri)
                last = stack.pop()
            else:
                break
        stack.append(last)
        stack.append(u)
    print(f'Stack after: {stack}')

print(f'\nFinal stack: {stack}')
print('Remaining:')
while len(stack) > 2:
    v1 = stack.pop()
    v2 = stack[-1]
    tri = (indices[sorted_idx[-1]], indices[v1], indices[v2])
    print(f'  Triangle: {tri}')
    tris.append(tri)

print(f'\nAll triangles: {tris}')
