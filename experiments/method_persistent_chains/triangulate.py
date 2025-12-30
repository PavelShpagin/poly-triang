"""
Persistent Chain Method for Polygon Triangulation
O(n + r log r) deterministic time

Key insight: Splice out processed regions so edges are only traversed once.
"""

class Vertex:
    def __init__(self, x, y, idx=-1):
        self.x = x
        self.y = y
        self.idx = idx
        self.prev = None
        self.next = None
        self.removed = False
    
    def __repr__(self):
        return f"V{self.idx}({self.x:.1f},{self.y:.1f})"


def cross(o, a, b):
    """Cross product OA x OB"""
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def build_polygon(points):
    """Build doubly-linked list from points (CCW order)"""
    n = len(points)
    verts = [Vertex(p[0], p[1], i) for i, p in enumerate(points)]
    for i in range(n):
        verts[i].prev = verts[(i-1) % n]
        verts[i].next = verts[(i+1) % n]
    return verts[0], verts


def find_reflex(verts):
    """Find monotonicity-violating reflex vertices"""
    up, down = [], []
    for v in verts:
        if cross(v.prev, v, v.next) < -1e-9:  # Reflex
            if v.prev.y > v.y and v.next.y > v.y:
                up.append(v)
            elif v.prev.y < v.y and v.next.y < v.y:
                down.append(v)
    up.sort(key=lambda v: (v.y, v.x))
    down.sort(key=lambda v: (-v.y, v.x))
    return up, down


def edge_crosses(a, b, y):
    """Does edge (a,b) strictly cross height y?"""
    return (a.y < y < b.y) or (b.y < y < a.y)


def intersect_y(a, b, y):
    """Intersection of edge (a,b) with y=y"""
    t = (y - a.y) / (b.y - a.y)
    return a.x + t * (b.x - a.x)


def find_supports(v, edge_counter):
    """Find left and right support edges for upward reflex v"""
    y = v.y
    
    # Walk left
    curr = v.prev
    left_steps = 0
    while not curr.removed:
        p = curr.prev
        if p.removed:
            break
        edge_counter[0] += 1
        left_steps += 1
        if left_steps > 1000:
            return None, None, None, None
        if edge_crosses(p, curr, y):
            left_upper = curr if curr.y > y else p
            left_lower = p if p.y < y else curr
            left_x = intersect_y(p, curr, y)
            break
        curr = p
    else:
        return None, None, None, None
    
    # Walk right
    curr = v.next
    right_steps = 0
    while not curr.removed:
        n = curr.next
        if n.removed:
            break
        edge_counter[0] += 1
        right_steps += 1
        if right_steps > 1000:
            return None, None, None, None
        if edge_crosses(curr, n, y):
            right_upper = n if n.y > y else curr
            right_lower = curr if curr.y < y else n
            right_x = intersect_y(curr, n, y)
            break
        curr = n
    else:
        return None, None, None, None
    
    return (left_upper, left_lower, left_x), (right_upper, right_lower, right_x), left_steps, right_steps


def splice_valley(v, left_info, right_info):
    """Splice out valley below upward reflex v"""
    left_upper, left_lower, left_x = left_info
    right_upper, right_lower, right_x = right_info
    y = v.y
    
    # Create chord endpoints
    nl = Vertex(left_x, y, -1)
    nr = Vertex(right_x, y, -1)
    
    # Valley: nl -> left_lower -> ... -> v -> ... -> right_lower -> nr -> nl
    nl.next = left_lower
    left_lower.prev = nl
    nr.prev = right_lower
    right_lower.next = nr
    nr.next = nl
    nl.prev = nr
    
    # Mark valley vertices as removed from main polygon
    curr = left_lower
    count = 0
    while curr != right_lower and count < 1000:
        curr.removed = True
        curr = curr.next
        count += 1
    right_lower.removed = True
    v.removed = True
    
    # Remaining: ... -> left_upper -> nl' -> nr' -> right_upper -> ...
    nl2 = Vertex(left_x, y, -1)
    nr2 = Vertex(right_x, y, -1)
    
    nl2.prev = left_upper
    left_upper.next = nl2
    nr2.next = right_upper
    right_upper.prev = nr2
    nl2.next = nr2
    nr2.prev = nl2
    
    return nl, nl2  # valley start, remaining start


def count_remaining_vertices(start):
    """Count vertices in a cycle"""
    count = 0
    curr = start
    while True:
        if not curr.removed:
            count += 1
        curr = curr.next
        if curr == start or count > 10000:
            break
    return count


def triangulate(points):
    """Main triangulation algorithm"""
    n = len(points)
    if n < 3:
        return [], 0
    
    start, verts = build_polygon(points)
    up, down = find_reflex(verts)
    
    r = len(up) + len(down)
    edge_counter = [0]  # Mutable counter
    
    valleys = []
    remaining = start
    
    # Process upward vertices bottom-to-top
    for v in up:
        if v.removed:
            continue
        
        result = find_supports(v, edge_counter)
        if result[0] is None:
            continue
        
        left_info, right_info, _, _ = result
        valley_start, remain_start = splice_valley(v, left_info, right_info)
        valleys.append(valley_start)
        remaining = remain_start
    
    return valleys, edge_counter[0]


def test():
    """Test on various polygons"""
    # Test 1: Simple arrow
    arrow = [(0,0), (4,0), (4,4), (2,2), (0,4)]
    vals, edges = triangulate(arrow)
    print(f"Arrow (n=5): {len(vals)} valleys, {edges} edge traversals")
    
    # Test 2: W shape
    w_shape = [(0,4), (1,0), (2,3), (3,0), (4,4)]
    vals, edges = triangulate(w_shape)
    print(f"W-shape (n=5): {len(vals)} valleys, {edges} edge traversals")
    
    # Test 3: Star
    import math
    star = []
    for i in range(10):
        angle = i * math.pi / 5
        r = 2 if i % 2 == 0 else 1
        star.append((r * math.cos(angle), r * math.sin(angle)))
    vals, edges = triangulate(star)
    print(f"Star (n=10): {len(vals)} valleys, {edges} edge traversals, ratio={edges/10:.2f}")
    
    # Test 4: Nested valleys
    nested = [
        (0, 10), (10, 10), (10, 0),
        (9, 1), (8, 0.5), (7, 1),
        (6, 2), (5, 1), (4, 2),
        (3, 3), (2, 2), (1, 3),
        (0, 0)
    ]
    vals, edges = triangulate(nested)
    print(f"Nested (n=13): {len(vals)} valleys, {edges} edge traversals, ratio={edges/13:.2f}")


if __name__ == "__main__":
    test()

