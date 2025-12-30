"""
Boundary Walk Triangulation with Horizon Tracking
Goal: O(n) deterministic time

Key idea: Process vertices in boundary order, maintain a "horizon" stack.
Clip ears when possible, defer when blocked by reflex vertices.
"""

from typing import List, Tuple, Set, Optional
from dataclasses import dataclass


@dataclass
class Point:
    x: float
    y: float
    idx: int = -1
    
    def __repr__(self):
        return f"P{self.idx}({self.x:.1f},{self.y:.1f})"


def cross(o: Point, a: Point, b: Point) -> float:
    """Cross product of vectors OA and OB"""
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def is_left_turn(a: Point, b: Point, c: Point) -> bool:
    """Check if a->b->c makes a left turn (CCW)"""
    return cross(a, b, c) > 1e-9


def is_reflex(prev: Point, curr: Point, next: Point) -> bool:
    """Check if curr is a reflex vertex (interior angle > 180)"""
    return cross(prev, curr, next) < -1e-9


def point_in_triangle(p: Point, a: Point, b: Point, c: Point) -> bool:
    """Check if point p is strictly inside triangle abc"""
    # Use barycentric coordinates
    d1 = cross(a, b, p)
    d2 = cross(b, c, p)
    d3 = cross(c, a, p)
    
    has_neg = (d1 < -1e-9) or (d2 < -1e-9) or (d3 < -1e-9)
    has_pos = (d1 > 1e-9) or (d2 > 1e-9) or (d3 > 1e-9)
    
    # Strictly inside if all same sign and non-zero
    if has_neg and has_pos:
        return False
    # Check not on boundary
    return abs(d1) > 1e-9 and abs(d2) > 1e-9 and abs(d3) > 1e-9


def mark_reflex_vertices(polygon: List[Point]) -> Set[int]:
    """Identify all reflex vertices"""
    n = len(polygon)
    reflex = set()
    for i in range(n):
        prev = polygon[(i - 1) % n]
        curr = polygon[i]
        next = polygon[(i + 1) % n]
        if is_reflex(prev, curr, next):
            reflex.add(i)
    return reflex


def can_clip(horizon: List[Point], v: Point, reflex_set: Set[int], 
             polygon: List[Point], processed: Set[int]) -> bool:
    """
    Check if we can clip the top of horizon to form a triangle with v.
    Returns True if triangle (w, u, v) is valid.
    """
    if len(horizon) < 2:
        return False
    
    u = horizon[-1]  # top
    w = horizon[-2]  # second from top
    
    # Must make a left turn (correct orientation for CCW polygon)
    if not is_left_turn(w, u, v):
        return False
    
    # Check that no unprocessed reflex vertex is inside the triangle
    for r_idx in reflex_set:
        if r_idx in processed:
            continue
        r = polygon[r_idx]
        # Skip vertices that are part of the triangle
        if r.idx == w.idx or r.idx == u.idx or r.idx == v.idx:
            continue
        if point_in_triangle(r, w, u, v):
            return False
    
    return True


def triangulate_boundary_walk(points: List[Tuple[float, float]]) -> List[Tuple[int, int, int]]:
    """
    Triangulate a simple polygon using boundary walk with horizon tracking.
    
    Args:
        points: List of (x, y) coordinates in CCW order
        
    Returns:
        List of triangles, each as (i, j, k) vertex indices
    """
    n = len(points)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]
    
    # Create Point objects with indices
    polygon = [Point(x, y, i) for i, (x, y) in enumerate(points)]
    
    # Phase 1: Mark reflex vertices
    reflex_set = mark_reflex_vertices(polygon)
    print(f"Polygon has {n} vertices, {len(reflex_set)} reflex")
    
    # Phase 2: Boundary walk with horizon tracking
    triangles = []
    horizon = []  # Stack of Points
    processed = set()  # Indices of fully processed vertices
    
    # Initialize with first two vertices
    horizon.append(polygon[0])
    horizon.append(polygon[1])
    
    # Process remaining vertices
    clip_attempts = 0
    successful_clips = 0
    
    for i in range(2, n):
        v = polygon[i]
        
        # Try to clip as many triangles as possible
        while can_clip(horizon, v, reflex_set, polygon, processed):
            clip_attempts += 1
            u = horizon.pop()
            w = horizon[-1]
            triangles.append((w.idx, u.idx, v.idx))
            successful_clips += 1
            processed.add(u.idx)
        
        clip_attempts += 1  # Count the failed attempt that exits the loop
        horizon.append(v)
    
    # Close the polygon: form fan from v_0 to remaining horizon
    while len(horizon) > 2:
        u = horizon.pop()
        w = horizon[-1]
        triangles.append((polygon[0].idx, w.idx, u.idx))
        successful_clips += 1
    
    print(f"Clip attempts: {clip_attempts}, successful: {successful_clips}")
    print(f"Triangles formed: {len(triangles)}")
    
    return triangles


def verify_triangulation(points: List[Tuple[float, float]], 
                         triangles: List[Tuple[int, int, int]]) -> bool:
    """Verify the triangulation is valid"""
    n = len(points)
    
    # Check triangle count
    if len(triangles) != n - 2:
        print(f"ERROR: Expected {n-2} triangles, got {len(triangles)}")
        return False
    
    # Check all triangles have correct orientation
    for t in triangles:
        i, j, k = t
        p1 = Point(points[i][0], points[i][1])
        p2 = Point(points[j][0], points[j][1])
        p3 = Point(points[k][0], points[k][1])
        if not is_left_turn(p1, p2, p3):
            print(f"ERROR: Triangle {t} has wrong orientation")
            return False
    
    print("Triangulation verified!")
    return True


def test_convex():
    """Test on convex polygon (should be trivial)"""
    print("\n=== Test: Convex Pentagon ===")
    import math
    pentagon = [(math.cos(2*math.pi*i/5), math.sin(2*math.pi*i/5)) for i in range(5)]
    triangles = triangulate_boundary_walk(pentagon)
    verify_triangulation(pentagon, triangles)
    print(f"Triangles: {triangles}")


def test_simple_nonconvex():
    """Test on simple non-convex polygon"""
    print("\n=== Test: Arrow Shape ===")
    arrow = [(0, 0), (4, 0), (4, 3), (2, 1.5), (0, 3)]
    triangles = triangulate_boundary_walk(arrow)
    verify_triangulation(arrow, triangles)
    print(f"Triangles: {triangles}")


def test_w_shape():
    """Test on W-shaped polygon with valleys"""
    print("\n=== Test: W Shape ===")
    w_shape = [(0, 3), (1, 0), (2, 2), (3, 0), (4, 3)]
    triangles = triangulate_boundary_walk(w_shape)
    verify_triangulation(w_shape, triangles)
    print(f"Triangles: {triangles}")


def test_nested_valleys():
    """Test with nested valleys - potential O(n^2) case"""
    print("\n=== Test: Nested Valleys ===")
    # Create polygon with multiple nested valleys
    nested = [
        (0, 10),   # 0: top-left
        (10, 10),  # 1: top-right
        (10, 0),   # 2: bottom-right
        (8, 2),    # 3: valley 1 right
        (6, 1),    # 4: valley 1 bottom (reflex)
        (4, 2),    # 5: valley 1 left
        (2, 3),    # 6: valley 2 right
        (1, 1),    # 7: valley 2 bottom (reflex)
        (0, 3),    # 8: valley 2 left
    ]
    # Close to form a simple polygon
    triangles = triangulate_boundary_walk(nested)
    verify_triangulation(nested, triangles)
    print(f"Triangles: {triangles}")


def test_star():
    """Test on star polygon"""
    print("\n=== Test: 5-Pointed Star ===")
    import math
    star = []
    for i in range(10):
        angle = i * math.pi / 5 - math.pi/2
        r = 2 if i % 2 == 0 else 1
        star.append((r * math.cos(angle), r * math.sin(angle)))
    triangles = triangulate_boundary_walk(star)
    verify_triangulation(star, triangles)
    print(f"Triangles: {triangles}")


def test_large_random():
    """Test on larger random polygon"""
    print("\n=== Test: Random Polygon (n=100) ===")
    import random
    random.seed(42)
    
    # Generate random points and compute convex hull ordering
    # (This ensures a simple polygon)
    n = 100
    angles = sorted([random.uniform(0, 2*3.14159) for _ in range(n)])
    
    # Create star-like polygon with random radii
    points = []
    for angle in angles:
        r = 1 + 0.5 * random.random()
        points.append((r * __import__('math').cos(angle), 
                      r * __import__('math').sin(angle)))
    
    triangles = triangulate_boundary_walk(points)
    verify_triangulation(points, triangles)
    
    # Count operations
    reflex_count = len(mark_reflex_vertices([Point(x, y, i) for i, (x, y) in enumerate(points)]))
    print(f"Reflex vertices: {reflex_count}")


def count_operations():
    """Count operations to verify O(n) behavior"""
    print("\n=== Operation Count Analysis ===")
    import math
    
    sizes = [10, 20, 50, 100, 200, 500]
    
    for n in sizes:
        # Create a polygon with many reflex vertices
        points = []
        for i in range(n):
            angle = 2 * math.pi * i / n
            # Alternate radii to create reflex vertices
            r = 2 if i % 3 == 0 else 1
            points.append((r * math.cos(angle), r * math.sin(angle)))
        
        polygon = [Point(x, y, i) for i, (x, y) in enumerate(points)]
        reflex_set = mark_reflex_vertices(polygon)
        r = len(reflex_set)
        
        # Count clip attempts
        horizon = [polygon[0], polygon[1]]
        processed = set()
        total_checks = 0
        reflex_checks = 0
        
        for i in range(2, n):
            v = polygon[i]
            while len(horizon) >= 2:
                total_checks += 1
                u = horizon[-1]
                w = horizon[-2]
                
                # Count reflex checks
                for r_idx in reflex_set:
                    if r_idx not in processed:
                        reflex_checks += 1
                
                if is_left_turn(w, u, v):
                    blocked = False
                    for r_idx in reflex_set:
                        if r_idx not in processed and r_idx not in {w.idx, u.idx, v.idx}:
                            if point_in_triangle(polygon[r_idx], w, u, v):
                                blocked = True
                                break
                    if not blocked:
                        horizon.pop()
                        processed.add(u.idx)
                        continue
                break
            horizon.append(v)
        
        print(f"n={n:4d}: reflex={r:3d}, total_checks={total_checks:5d}, "
              f"reflex_checks={reflex_checks:6d}, ratio={reflex_checks/n:.2f}")


if __name__ == "__main__":
    test_convex()
    test_simple_nonconvex()
    test_w_shape()
    test_nested_valleys()
    test_star()
    test_large_random()
    count_operations()

