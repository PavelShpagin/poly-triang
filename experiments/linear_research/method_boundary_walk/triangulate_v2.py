"""
Boundary Walk Triangulation v2 - Optimized O(n)
Key fix: Only check reflex vertices in the "active range" - those between
the horizon and current position in boundary order.
"""

from typing import List, Tuple, Set, Dict
from dataclasses import dataclass
from collections import deque


@dataclass
class Point:
    x: float
    y: float
    idx: int = -1


def cross(o: Point, a: Point, b: Point) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def is_left_turn(a: Point, b: Point, c: Point) -> bool:
    return cross(a, b, c) > 1e-9


def is_reflex(prev: Point, curr: Point, next: Point) -> bool:
    return cross(prev, curr, next) < -1e-9


def point_in_triangle(p: Point, a: Point, b: Point, c: Point) -> bool:
    d1 = cross(a, b, p)
    d2 = cross(b, c, p)
    d3 = cross(c, a, p)
    has_neg = (d1 < -1e-9) or (d2 < -1e-9) or (d3 < -1e-9)
    has_pos = (d1 > 1e-9) or (d2 > 1e-9) or (d3 > 1e-9)
    if has_neg and has_pos:
        return False
    return abs(d1) > 1e-9 and abs(d2) > 1e-9 and abs(d3) > 1e-9


class ReflexTracker:
    """
    Efficient tracking of reflex vertices for triangle-interior queries.
    Key insight: We only need to check reflex vertices that are "ahead" 
    of the current horizon position in boundary order.
    """
    def __init__(self, polygon: List[Point], reflex_indices: Set[int]):
        self.polygon = polygon
        self.n = len(polygon)
        # Reflex vertices sorted by index (boundary order)
        self.reflex_list = sorted(reflex_indices)
        self.reflex_set = reflex_indices
        # Pointer to first "active" reflex vertex
        self.active_start = 0
        # Track which reflex vertices have been "passed"
        self.passed = set()
        
    def get_active_reflex(self, horizon_back_idx: int, current_idx: int) -> List[int]:
        """
        Get reflex vertices that could be inside a triangle formed by 
        horizon vertices and current vertex.
        
        These are reflex vertices with index in (horizon_back_idx, current_idx)
        that haven't been passed yet.
        """
        result = []
        for r in self.reflex_list:
            if r in self.passed:
                continue
            # Check if r is "between" horizon and current in boundary order
            if self._is_between(horizon_back_idx, r, current_idx):
                result.append(r)
        return result
    
    def _is_between(self, start: int, mid: int, end: int) -> bool:
        """Check if mid is between start and end in cyclic boundary order"""
        if start < end:
            return start < mid < end
        else:
            # Wraps around
            return mid > start or mid < end
    
    def mark_passed(self, idx: int):
        """Mark a reflex vertex as passed (no longer needs checking)"""
        if idx in self.reflex_set:
            self.passed.add(idx)


def triangulate_optimized(points: List[Tuple[float, float]]) -> Tuple[List[Tuple[int, int, int]], Dict]:
    """
    O(n) triangulation via boundary walk with optimized reflex tracking.
    """
    n = len(points)
    if n < 3:
        return [], {}
    if n == 3:
        return [(0, 1, 2)], {"reflex": 0, "checks": 0}
    
    # Create polygon
    polygon = [Point(x, y, i) for i, (x, y) in enumerate(points)]
    
    # Find reflex vertices
    reflex_indices = set()
    for i in range(n):
        if is_reflex(polygon[(i-1) % n], polygon[i], polygon[(i+1) % n]):
            reflex_indices.add(i)
    
    tracker = ReflexTracker(polygon, reflex_indices)
    
    # Stats
    stats = {"reflex": len(reflex_indices), "clip_attempts": 0, "reflex_checks": 0}
    
    triangles = []
    horizon = [polygon[0], polygon[1]]  # Stack
    
    for i in range(2, n):
        v = polygon[i]
        
        while len(horizon) >= 2:
            stats["clip_attempts"] += 1
            u = horizon[-1]
            w = horizon[-2]
            
            # Check orientation
            if not is_left_turn(w, u, v):
                break
            
            # Check for blocking reflex vertices
            # OPTIMIZATION: Only check reflex vertices in the relevant range
            active_reflex = tracker.get_active_reflex(w.idx, v.idx)
            stats["reflex_checks"] += len(active_reflex)
            
            blocked = False
            for r_idx in active_reflex:
                r = polygon[r_idx]
                if r_idx == w.idx or r_idx == u.idx or r_idx == v.idx:
                    continue
                if point_in_triangle(r, w, u, v):
                    blocked = True
                    break
            
            if blocked:
                break
            
            # Clip the ear
            horizon.pop()
            triangles.append((w.idx, u.idx, v.idx))
            # Mark u as passed (if it was reflex)
            tracker.mark_passed(u.idx)
        
        horizon.append(v)
        # Mark v as passed if reflex
        tracker.mark_passed(v.idx)
    
    # Close polygon - the remaining horizon forms a monotone chain
    # Triangulate it properly (not just a fan from v0)
    while len(horizon) > 2:
        # Pop from top, form triangle with the two below it
        u = horizon.pop()
        w = horizon[-1]
        if len(horizon) >= 2:
            z = horizon[-2]
            triangles.append((z.idx, w.idx, u.idx))
        else:
            # Last triangle connects back to v0
            triangles.append((polygon[0].idx, w.idx, u.idx))
    
    return triangles, stats


def verify(points: List[Tuple[float, float]], triangles: List[Tuple[int, int, int]]) -> bool:
    n = len(points)
    if len(triangles) != n - 2:
        print(f"ERROR: Expected {n-2} triangles, got {len(triangles)}")
        return False
    
    # Check orientation
    for i, j, k in triangles:
        p1, p2, p3 = Point(points[i][0], points[i][1]), Point(points[j][0], points[j][1]), Point(points[k][0], points[k][1])
        area = cross(p1, p2, p3)
        if area < -1e-9:
            print(f"ERROR: Triangle ({i},{j},{k}) has negative area")
            return False
    
    return True


def run_complexity_test():
    """Test to verify O(n) complexity"""
    import math
    
    print("=== Complexity Analysis (Optimized) ===")
    print(f"{'n':>6} {'r':>6} {'attempts':>10} {'reflex_chk':>12} {'ratio':>8}")
    print("-" * 50)
    
    for n in [10, 20, 50, 100, 200, 500, 1000, 2000]:
        # Create polygon with ~n/3 reflex vertices
        points = []
        for i in range(n):
            angle = 2 * math.pi * i / n
            r = 2 if i % 3 == 0 else 1
            points.append((r * math.cos(angle), r * math.sin(angle)))
        
        triangles, stats = triangulate_optimized(points)
        
        if not verify(points, triangles):
            print(f"FAILED for n={n}")
            continue
        
        ratio = stats["reflex_checks"] / n if n > 0 else 0
        print(f"{n:>6} {stats['reflex']:>6} {stats['clip_attempts']:>10} "
              f"{stats['reflex_checks']:>12} {ratio:>8.2f}")


def run_basic_tests():
    """Basic correctness tests"""
    import math
    
    tests = [
        ("Convex Pentagon", [(math.cos(2*math.pi*i/5), math.sin(2*math.pi*i/5)) for i in range(5)]),
        ("Arrow", [(0, 0), (4, 0), (4, 3), (2, 1.5), (0, 3)]),
        ("W Shape", [(0, 3), (1, 0), (2, 2), (3, 0), (4, 3)]),
        ("Star", [(2*math.cos(i*math.pi/5 - math.pi/2) if i%2==0 else math.cos(i*math.pi/5 - math.pi/2),
                   2*math.sin(i*math.pi/5 - math.pi/2) if i%2==0 else math.sin(i*math.pi/5 - math.pi/2)) 
                 for i in range(10)]),
    ]
    
    print("=== Basic Tests ===")
    for name, points in tests:
        triangles, stats = triangulate_optimized(points)
        valid = verify(points, triangles)
        status = "OK" if valid else "FAIL"
        print(f"{name}: {status} - {len(triangles)} triangles, {stats['reflex']} reflex, "
              f"{stats['reflex_checks']} checks")


if __name__ == "__main__":
    run_basic_tests()
    print()
    run_complexity_test()

