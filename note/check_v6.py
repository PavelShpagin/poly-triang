#!/usr/bin/env python3
pts = [
    (1.5, 1.5),  # v0
    (3.0, 0.0),  # v1
    (5.0, 2.5),  # v2
    (8.0, 1.5),  # v3
    (6.5, 3.5),  # v4
    (8.0, 5.5),  # v5
    (7.0, 7.0),  # v6
    (5.5, 5.0),  # v7
    (4.0, 6.5),  # v8
    (2.5, 4.0),  # v9
    (1.0, 5.5),  # v10
    (0.0, 2.5),  # v11
]
n = len(pts)

def prv(i): return (i - 1 + n) % n
def nxt(i): return (i + 1) % n

# v6 at (7.0, 7.0)
v = 6
p, nx = prv(v), nxt(v)
print(f"v{v} at {pts[v]}")
print(f"  prv({v}) = v{p} at {pts[p]}")
print(f"  nxt({v}) = v{nx} at {pts[nx]}")
print(f"  prv is below: {pts[p][1] < pts[v][1]}")
print(f"  nxt is below: {pts[nx][1] < pts[v][1]}")
print()

# For END vertex, both neighbors should be ABOVE (greater y)
# But v6 is at y=7.0, prv(6)=v5 at y=5.5, nxt(6)=v7 at y=5.0
# Both are BELOW v6, so v6 should be START, not END!

print("Classification check:")
for i in range(n):
    p, nx = prv(i), nxt(i)
    pb = pts[p][1] < pts[i][1]  # prev below
    nb = pts[nx][1] < pts[i][1]  # next below
    
    # Cross product for reflex check
    ax, ay = pts[p]
    bx, by = pts[i]
    cx, cy = pts[nx]
    cross = (bx - ax) * (cy - by) - (by - ay) * (cx - bx)
    
    if not pb and not nb:  # both above -> local min
        vtype = "START" if cross > 0 else "SPLIT"
    elif pb and nb:  # both below -> local max
        vtype = "END" if cross > 0 else "MERGE"
    else:
        vtype = "REG_R" if pb else "REG_L"
    
    print(f"v{i} at {pts[i]}: pb={pb}, nb={nb}, cross={cross:.2f} -> {vtype}")

