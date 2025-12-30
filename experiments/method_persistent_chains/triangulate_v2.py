"""
Persistent Chain Method v2 - Correct Implementation
Key fix: Actually modify the linked list structure after each splice
"""

class Vertex:
    __slots__ = ['x', 'y', 'idx', 'prev', 'next']
    def __init__(self, x, y, idx=-1):
        self.x, self.y, self.idx = x, y, idx
        self.prev = self.next = None


def cross(o, a, b):
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def build_polygon(points):
    n = len(points)
    verts = [Vertex(p[0], p[1], i) for i, p in enumerate(points)]
    for i in range(n):
        verts[i].prev = verts[(i-1) % n]
        verts[i].next = verts[(i+1) % n]
    return verts[0], verts


def find_reflex(start):
    """Find reflex vertices from current polygon state"""
    up, down = [], []
    curr = start
    seen = set()
    while id(curr) not in seen:
        seen.add(id(curr))
        if cross(curr.prev, curr, curr.next) < -1e-9:
            if curr.prev.y > curr.y and curr.next.y > curr.y:
                up.append(curr)
            elif curr.prev.y < curr.y and curr.next.y < curr.y:
                down.append(curr)
        curr = curr.next
        if len(seen) > 10000:
            break
    up.sort(key=lambda v: (v.y, v.x))
    down.sort(key=lambda v: (-v.y, v.x))
    return up, down


def edge_crosses(a, b, y):
    return (a.y < y < b.y) or (b.y < y < a.y)


def intersect_y(a, b, y):
    if abs(b.y - a.y) < 1e-12:
        return (a.x + b.x) / 2
    t = (y - a.y) / (b.y - a.y)
    return a.x + t * (b.x - a.x)


def find_left_support(v, counter, max_steps=500):
    """Find left support edge starting from v.prev"""
    y = v.y
    curr = v.prev
    steps = 0
    while steps < max_steps:
        p = curr.prev
        counter[0] += 1
        steps += 1
        if edge_crosses(p, curr, y):
            if curr.y > y:
                return curr, p, intersect_y(p, curr, y)
            else:
                return p, curr, intersect_y(curr, p, y)
        curr = p
        if curr == v:
            break
    return None, None, None


def find_right_support(v, counter, max_steps=500):
    """Find right support edge starting from v.next"""
    y = v.y
    curr = v.next
    steps = 0
    while steps < max_steps:
        n = curr.next
        counter[0] += 1
        steps += 1
        if edge_crosses(curr, n, y):
            if curr.y > y:
                return curr, n, intersect_y(curr, n, y)
            else:
                return n, curr, intersect_y(n, curr, y)
        curr = n
        if curr == v:
            break
    return None, None, None


def splice_upward(v, lu, ll, lx, ru, rl, rx):
    """
    Splice polygon at upward reflex v.
    
    lu, ll = left upper, left lower vertices
    ru, rl = right upper, right lower vertices
    lx, rx = x-coordinates of chord endpoints
    
    Returns: (valley_start, remaining_start)
    """
    y = v.y
    
    # Create chord endpoint vertices
    pl = Vertex(lx, y)  # left point
    pr = Vertex(rx, y)  # right point
    
    # Valley cycle: pl -> ll -> ... -> v -> ... -> rl -> pr -> pl
    # First, find the chain from ll to rl through v
    # This is the "below" part
    
    # Connect pl into valley
    pl.next = ll
    ll.prev = pl
    
    # Connect pr into valley
    pr.prev = rl
    rl.next = pr
    
    # Close valley
    pr.next = pl
    pl.prev = pr
    
    # Remaining polygon: ... -> lu -> pl' -> pr' -> ru -> ...
    pl2 = Vertex(lx, y)
    pr2 = Vertex(rx, y)
    
    # Connect into remaining
    pl2.prev = lu
    lu.next = pl2
    
    pr2.next = ru
    ru.prev = pr2
    
    # Connect chord in remaining
    pl2.next = pr2
    pr2.prev = pl2
    
    return pl, pl2


def count_vertices(start, limit=1000):
    """Count vertices in cycle"""
    count = 0
    curr = start
    while count < limit:
        count += 1
        curr = curr.next
        if curr == start:
            break
    return count


def process_upward(start, counter):
    """Process all upward reflex vertices"""
    valleys = []
    remaining = start
    
    # Re-find reflex vertices from current state
    up, _ = find_reflex(remaining)
    
    for v in up:
        # Check v is still in remaining polygon
        curr = remaining
        found = False
        for _ in range(count_vertices(remaining)):
            if curr == v:
                found = True
                break
            curr = curr.next
        if not found:
            continue
        
        # Find supports
        lu, ll, lx = find_left_support(v, counter)
        if lu is None:
            continue
        
        ru, rl, rx = find_right_support(v, counter)
        if ru is None:
            continue
        
        # Splice
        valley, remain = splice_upward(v, lu, ll, lx, ru, rl, rx)
        valleys.append(valley)
        remaining = remain
    
    return valleys, remaining


def triangulate(points):
    """Main triangulation"""
    n = len(points)
    if n < 3:
        return [], 0
    
    start, verts = build_polygon(points)
    counter = [0]
    
    valleys, remaining = process_upward(start, counter)
    
    return valleys, counter[0]


def test():
    # Test 1: Simple
    simple = [(0,0), (2,0), (2,2), (1,1), (0,2)]
    vals, edges = triangulate(simple)
    print(f"Simple (n=5): {len(vals)} valleys, {edges} edges, ratio={edges/5:.2f}")
    
    # Test 2: Nested valleys - this is the critical test
    nested = []
    for i in range(5):
        # Create V shapes at increasing heights
        nested.extend([
            (i, 5-i),      # left wall
            (i+0.5, 4-i),  # valley bottom
            (i+1, 5-i),    # right wall
        ])
    # Close the polygon
    nested.append((5, 5))
    nested.append((0, 5))
    
    vals, edges = triangulate(nested)
    n = len(nested)
    print(f"Nested (n={n}): {len(vals)} valleys, {edges} edges, ratio={edges/n:.2f}")
    
    # Test 3: Simple nested
    simple_nested = [
        (0, 4), (4, 4),  # top
        (4, 0),          # right
        (3, 1), (2, 0.5), (1, 1),  # inner valley
        (0, 0),          # left
    ]
    vals, edges = triangulate(simple_nested)
    n = len(simple_nested)
    print(f"Simple nested (n={n}): {len(vals)} valleys, {edges} edges, ratio={edges/n:.2f}")


if __name__ == "__main__":
    test()

