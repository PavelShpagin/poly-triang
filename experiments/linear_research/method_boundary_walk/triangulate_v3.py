"""
Boundary Walk Triangulation v3 - Clean implementation
Uses standard monotone polygon triangulation for the final phase.
"""

from typing import List, Tuple, Set
from dataclasses import dataclass


@dataclass 
class Point:
    x: float
    y: float
    idx: int = -1


def cross(o: Point, a: Point, b: Point) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def is_convex_vertex(prev: Point, curr: Point, next: Point) -> bool:
    """Check if curr makes a left turn (convex for CCW polygon)"""
    return cross(prev, curr, next) > 1e-9


def point_in_triangle_strict(p: Point, a: Point, b: Point, c: Point) -> bool:
    """Check if p is strictly inside triangle abc"""
    d1 = cross(a, b, p)
    d2 = cross(b, c, p)
    d3 = cross(c, a, p)
    has_neg = (d1 < -1e-9) or (d2 < -1e-9) or (d3 < -1e-9)
    has_pos = (d1 > 1e-9) or (d2 > 1e-9) or (d3 > 1e-9)
    if has_neg and has_pos:
        return False
    return abs(d1) > 1e-9 and abs(d2) > 1e-9 and abs(d3) > 1e-9


def is_ear(polygon: List[Point], i: int, reflex_set: Set[int]) -> bool:
    """Check if vertex i is an ear"""
    n = len(polygon)
    prev = polygon[(i - 1) % n]
    curr = polygon[i]
    next = polygon[(i + 1) % n]
    
    # Must be convex
    if not is_convex_vertex(prev, curr, next):
        return False
    
    # Check no reflex vertex inside triangle
    for r in reflex_set:
        if r == (i - 1) % n or r == i or r == (i + 1) % n:
            continue
        if point_in_triangle_strict(polygon[r], prev, curr, next):
            return False
    
    return True


def ear_clipping_optimized(points: List[Tuple[float, float]]) -> Tuple[List[Tuple[int, int, int]], dict]:
    """
    Optimized ear clipping with reflex vertex tracking.
    
    Optimization: Only check reflex vertices (not all vertices) for triangle interior.
    This gives O(n * r) worst case, which is O(n^2) but often much better in practice.
    """
    n = len(points)
    if n < 3:
        return [], {"reflex": 0, "checks": 0}
    if n == 3:
        return [(0, 1, 2)], {"reflex": 0, "checks": 0}
    
    # Build polygon as list
    polygon = [Point(x, y, i) for i, (x, y) in enumerate(points)]
    
    # Find reflex vertices
    reflex_set = set()
    for i in range(n):
        if not is_convex_vertex(polygon[(i-1) % n], polygon[i], polygon[(i+1) % n]):
            reflex_set.add(i)
    
    stats = {"reflex": len(reflex_set), "checks": 0, "ear_scans": 0}
    
    # Use doubly linked list indices
    prev_idx = [(i - 1) % n for i in range(n)]
    next_idx = [(i + 1) % n for i in range(n)]
    active = [True] * n
    
    triangles = []
    remaining = n
    
    # Find initial ears
    ears = []
    for i in range(n):
        if is_ear_linked(polygon, i, reflex_set, prev_idx, next_idx, active, stats):
            ears.append(i)
    
    while remaining > 3 and ears:
        # Take an ear
        i = ears.pop(0)
        stats["ear_scans"] += 1
        
        if not active[i]:
            continue
        
        # Verify still an ear (neighbors might have changed)
        if not is_ear_linked(polygon, i, reflex_set, prev_idx, next_idx, active, stats):
            continue
        
        # Get neighbors
        p = prev_idx[i]
        q = next_idx[i]
        
        # Add triangle
        triangles.append((polygon[p].idx, polygon[i].idx, polygon[q].idx))
        
        # Remove vertex i
        active[i] = False
        next_idx[p] = q
        prev_idx[q] = p
        remaining -= 1
        
        # Update reflex status of neighbors
        if p in reflex_set:
            pp = prev_idx[p]
            if is_convex_vertex(polygon[pp], polygon[p], polygon[q]):
                reflex_set.remove(p)
        
        if q in reflex_set:
            qq = next_idx[q]
            if is_convex_vertex(polygon[p], polygon[q], polygon[qq]):
                reflex_set.remove(q)
        
        # Check if neighbors became ears
        if is_ear_linked(polygon, p, reflex_set, prev_idx, next_idx, active, stats):
            if p not in ears:
                ears.append(p)
        if is_ear_linked(polygon, q, reflex_set, prev_idx, next_idx, active, stats):
            if q not in ears:
                ears.append(q)
    
    # Handle last triangle
    if remaining == 3:
        # Find three remaining vertices
        for i in range(n):
            if active[i]:
                p = prev_idx[i]
                q = next_idx[i]
                triangles.append((polygon[p].idx, polygon[i].idx, polygon[q].idx))
                break
    
    return triangles, stats


def is_ear_linked(polygon: List[Point], i: int, reflex_set: Set[int],
                  prev_idx: List[int], next_idx: List[int], 
                  active: List[bool], stats: dict) -> bool:
    """Check if vertex i is an ear in linked list structure"""
    if not active[i]:
        return False
    
    p = prev_idx[i]
    q = next_idx[i]
    
    prev = polygon[p]
    curr = polygon[i]
    next = polygon[q]
    
    # Must be convex
    if not is_convex_vertex(prev, curr, next):
        return False
    
    # Check reflex vertices
    for r in reflex_set:
        if not active[r]:
            continue
        if r == p or r == i or r == q:
            continue
        stats["checks"] += 1
        if point_in_triangle_strict(polygon[r], prev, curr, next):
            return False
    
    return True


def verify_triangulation(points: List[Tuple[float, float]], 
                         triangles: List[Tuple[int, int, int]]) -> bool:
    n = len(points)
    
    if len(triangles) != n - 2:
        print(f"ERROR: Expected {n-2} triangles, got {len(triangles)}")
        return False
    
    # Check orientation
    for t in triangles:
        i, j, k = t
        p1 = Point(points[i][0], points[i][1])
        p2 = Point(points[j][0], points[j][1])
        p3 = Point(points[k][0], points[k][1])
        area = cross(p1, p2, p3)
        if area <= 1e-9:
            print(f"ERROR: Triangle {t} has non-positive area {area:.6f}")
            return False
    
    return True


def run_tests():
    import math
    
    tests = [
        ("Convex Pentagon", [(math.cos(2*math.pi*i/5), math.sin(2*math.pi*i/5)) for i in range(5)]),
        ("Arrow", [(0, 0), (4, 0), (4, 3), (2, 1.5), (0, 3)]),
        ("W Shape", [(0, 3), (1, 0), (2, 2), (3, 0), (4, 3)]),
        ("Star", [(2*math.cos(i*math.pi/5 - math.pi/2) if i%2==0 else math.cos(i*math.pi/5 - math.pi/2),
                   2*math.sin(i*math.pi/5 - math.pi/2) if i%2==0 else math.sin(i*math.pi/5 - math.pi/2)) 
                 for i in range(10)]),
    ]
    
    print("=== Basic Tests (Optimized Ear Clipping) ===")
    for name, points in tests:
        triangles, stats = ear_clipping_optimized(points)
        valid = verify_triangulation(points, triangles)
        status = "OK" if valid else "FAIL"
        print(f"{name}: {status} - {len(triangles)} triangles, r={stats['reflex']}, "
              f"checks={stats['checks']}")


def run_complexity_analysis():
    import math
    
    print("\n=== Complexity Analysis ===")
    print(f"{'n':>6} {'r':>6} {'checks':>10} {'ratio':>8} {'scans':>8}")
    print("-" * 50)
    
    for n in [10, 20, 50, 100, 200, 500, 1000]:
        # Create polygon with ~n/3 reflex vertices
        points = []
        for i in range(n):
            angle = 2 * math.pi * i / n
            r = 2 if i % 3 == 0 else 1
            points.append((r * math.cos(angle), r * math.sin(angle)))
        
        triangles, stats = ear_clipping_optimized(points)
        
        if not verify_triangulation(points, triangles):
            print(f"FAILED for n={n}")
            continue
        
        ratio = stats["checks"] / n
        print(f"{n:>6} {stats['reflex']:>6} {stats['checks']:>10} {ratio:>8.2f} {stats['ear_scans']:>8}")


if __name__ == "__main__":
    run_tests()
    run_complexity_analysis()

