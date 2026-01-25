import math

def orient(ax, ay, bx, by, cx, cy) -> float:
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)

def on_segment(ax, ay, bx, by, cx, cy) -> bool:
    return (
        min(ax, bx) - 1e-10 <= cx <= max(ax, bx) + 1e-10
        and min(ay, by) - 1e-10 <= cy <= max(ay, by) + 1e-10
    )

EPS = 1e-12
def segments_properly_intersect(a, b, c, d) -> bool:
    ax, ay = a
    bx, by = b
    cx, cy = c
    dx, dy = d

    o1 = orient(ax, ay, bx, by, cx, cy)
    o2 = orient(ax, ay, bx, by, dx, dy)
    o3 = orient(cx, cy, dx, dy, ax, ay)
    o4 = orient(cx, cy, dx, dy, bx, by)

    def sgn(x):
        if x > EPS: return 1
        if x < -EPS: return -1
        return 0

    s1, s2, s3, s4 = sgn(o1), sgn(o2), sgn(o3), sgn(o4)

    if s1 * s2 < 0 and s3 * s4 < 0:
        return True

    if s1 == 0 and on_segment(ax, ay, bx, by, cx, cy): return True
    if s2 == 0 and on_segment(ax, ay, bx, by, dx, dy): return True
    if s3 == 0 and on_segment(cx, cy, dx, dy, ax, ay): return True
    if s4 == 0 and on_segment(cx, cy, dx, dy, bx, by): return True

    return False

pts = []
with open("paper/CGAT/CGATv2/fail.poly") as f:
    lines = f.readlines()
    for line in lines[1:]: # skip count
        pts.append(tuple(map(float, line.split())))

p239 = pts[239]
p240 = pts[240]
p241 = pts[241]
p242 = pts[242]

print(f"239: {p239}")
print(f"240: {p240}")
print(f"241: {p241}")
print(f"242: {p242}")

print(f"Intersect (239,242) vs (240,241): {segments_properly_intersect(p239, p242, p240, p241)}")
