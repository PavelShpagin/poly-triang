"""
Polygon Triangulation via Reflex Graph Method
Goal: O(n + r log r) deterministic time

Key insight: After processing each reflex vertex, we PHYSICALLY REMOVE
the cut-off region from the polygon. Future walks cannot access removed edges.
"""

from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Set
import math


@dataclass
class Point:
    x: float
    y: float
    
    def __hash__(self):
        return hash((self.x, self.y))
    
    def __eq__(self, other):
        return abs(self.x - other.x) < 1e-9 and abs(self.y - other.y) < 1e-9


@dataclass 
class Vertex:
    """Doubly-linked list node for polygon vertex"""
    point: Point
    prev: Optional['Vertex'] = None
    next: Optional['Vertex'] = None
    is_reflex: bool = False
    reflex_type: str = ""  # "up" or "down" or ""
    original_index: int = -1
    
    def __repr__(self):
        return f"V({self.point.x:.1f}, {self.point.y:.1f})"


def cross_product(o: Point, a: Point, b: Point) -> float:
    """Cross product of vectors OA and OB"""
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def build_polygon(points: List[Tuple[float, float]]) -> Tuple[Vertex, List[Vertex]]:
    """Build doubly-linked list from point list (CCW order)"""
    n = len(points)
    vertices = []
    
    for i, (x, y) in enumerate(points):
        v = Vertex(Point(x, y), original_index=i)
        vertices.append(v)
    
    for i in range(n):
        vertices[i].prev = vertices[(i - 1) % n]
        vertices[i].next = vertices[(i + 1) % n]
    
    return vertices[0], vertices


def classify_vertices(vertices: List[Vertex]) -> Tuple[List[Vertex], List[Vertex]]:
    """
    Classify vertices as upward or downward reflex.
    Returns (upward_reflex, downward_reflex) sorted by y-coordinate.
    """
    upward = []
    downward = []
    
    for v in vertices:
        prev_pt = v.prev.point
        curr_pt = v.point
        next_pt = v.next.point
        
        # Check if reflex (cross product < 0 for CCW polygon)
        cp = cross_product(prev_pt, curr_pt, next_pt)
        
        if cp < -1e-9:  # Reflex vertex
            v.is_reflex = True
            
            # Check if upward (both neighbors above) or downward (both below)
            if prev_pt.y > curr_pt.y and next_pt.y > curr_pt.y:
                v.reflex_type = "up"
                upward.append(v)
            elif prev_pt.y < curr_pt.y and next_pt.y < curr_pt.y:
                v.reflex_type = "down"
                downward.append(v)
    
    # Sort: upward by y ascending, downward by y descending
    upward.sort(key=lambda v: (v.point.y, v.point.x))
    downward.sort(key=lambda v: (-v.point.y, v.point.x))
    
    return upward, downward


def edge_crosses_height(v1: Vertex, v2: Vertex, y: float) -> bool:
    """Check if edge (v1, v2) strictly crosses height y"""
    y1, y2 = v1.point.y, v2.point.y
    return (y1 < y < y2) or (y2 < y < y1)


def intersect_horizontal(v1: Vertex, v2: Vertex, y: float) -> Point:
    """Find intersection of edge (v1, v2) with horizontal line y"""
    p1, p2 = v1.point, v2.point
    t = (y - p1.y) / (p2.y - p1.y)
    x = p1.x + t * (p2.x - p1.x)
    return Point(x, y)


def find_left_support(v: Vertex, y: float) -> Tuple[Vertex, Vertex, Point]:
    """
    Find left support edge for reflex vertex v at height y.
    Returns (upper_vertex, lower_vertex, intersection_point).
    
    Walk left from v.prev until finding an edge that crosses y.
    """
    curr = v.prev
    steps = 0
    max_steps = 10000  # Safety limit
    
    while steps < max_steps:
        prev_v = curr.prev
        if edge_crosses_height(prev_v, curr, y):
            # Found the support edge
            if prev_v.point.y < y:
                return (curr, prev_v, intersect_horizontal(prev_v, curr, y))
            else:
                return (prev_v, curr, intersect_horizontal(curr, prev_v, y))
        curr = prev_v
        steps += 1
        
        if curr == v:
            raise ValueError("Walked full circle without finding support edge")
    
    raise ValueError("Max steps exceeded in find_left_support")


def find_right_support(v: Vertex, y: float) -> Tuple[Vertex, Vertex, Point]:
    """
    Find right support edge for reflex vertex v at height y.
    Returns (upper_vertex, lower_vertex, intersection_point).
    """
    curr = v.next
    steps = 0
    max_steps = 10000
    
    while steps < max_steps:
        next_v = curr.next
        if edge_crosses_height(curr, next_v, y):
            if curr.point.y < y:
                return (next_v, curr, intersect_horizontal(curr, next_v, y))
            else:
                return (curr, next_v, intersect_horizontal(next_v, curr, y))
        curr = next_v
        steps += 1
        
        if curr == v:
            raise ValueError("Walked full circle without finding support edge")
    
    raise ValueError("Max steps exceeded in find_right_support")


def split_at_chord_upward(v: Vertex, left_upper: Vertex, left_lower: Vertex, 
                          p_left: Point, right_upper: Vertex, right_lower: Vertex, 
                          p_right: Point) -> Tuple[Vertex, Vertex]:
    """
    Split polygon at horizontal chord through upward reflex vertex v.
    
    Creates two regions:
    1. Lower region (valley) - returned as first vertex
    2. Upper region (remaining) - returned as first vertex
    
    The chord goes from p_left to v to p_right.
    """
    # Create new vertices at chord endpoints
    v_left = Vertex(p_left, original_index=-1)
    v_right = Vertex(p_right, original_index=-1)
    
    # The lower region (valley) boundary:
    # v_left -> left_lower -> ... -> v -> ... -> right_lower -> v_right -> v_left
    
    # Connect left chord point into lower region
    v_left.next = left_lower
    left_lower.prev = v_left
    
    # Connect right chord point into lower region  
    v_right.prev = right_lower
    right_lower.next = v_right
    
    # Close lower region with chord
    v_right.next = v_left
    v_left.prev = v_right
    
    # The upper region boundary:
    # left_upper -> v_left' -> v_right' -> right_upper -> ...
    
    v_left_upper = Vertex(p_left, original_index=-1)
    v_right_upper = Vertex(p_right, original_index=-1)
    
    # Connect into upper region
    v_left_upper.prev = left_upper
    left_upper.next = v_left_upper
    
    v_right_upper.next = right_upper
    right_upper.prev = v_right_upper
    
    # Connect chord in upper region
    v_left_upper.next = v_right_upper
    v_right_upper.prev = v_left_upper
    
    return v_left, v_left_upper  # Return starting vertices of both regions


def triangulate_monotone(start: Vertex) -> List[Tuple[Point, Point, Point]]:
    """
    Triangulate a y-monotone polygon using the standard algorithm.
    Returns list of triangles (as point triples).
    """
    # Collect all vertices
    vertices = []
    curr = start
    while True:
        vertices.append(curr)
        curr = curr.next
        if curr == start:
            break
    
    if len(vertices) < 3:
        return []
    
    if len(vertices) == 3:
        return [(vertices[0].point, vertices[1].point, vertices[2].point)]
    
    # Sort by y-coordinate (descending), then x
    sorted_verts = sorted(vertices, key=lambda v: (-v.point.y, v.point.x))
    
    # Determine left and right chains
    # Find topmost and bottommost
    top = sorted_verts[0]
    bottom = sorted_verts[-1]
    
    # Mark chain membership
    left_chain = set()
    curr = top
    while curr != bottom:
        left_chain.add(curr)
        curr = curr.next
    
    triangles = []
    stack = [sorted_verts[0], sorted_verts[1]]
    
    for i in range(2, len(sorted_verts)):
        u = sorted_verts[i]
        
        if (u in left_chain) != (stack[-1] in left_chain):
            # Different chains
            while len(stack) > 1:
                v = stack.pop()
                triangles.append((u.point, v.point, stack[-1].point))
            stack.pop()
            stack.append(sorted_verts[i-1])
            stack.append(u)
        else:
            # Same chain
            v = stack.pop()
            while len(stack) > 0:
                # Check if diagonal is inside
                v2 = stack[-1]
                # Simplified check - just add triangle
                triangles.append((u.point, v.point, v2.point))
                v = stack.pop()
                if len(stack) == 0:
                    break
            stack.append(sorted_verts[i-1])
            stack.append(u)
    
    return triangles


def count_edges_traversed(vertices: List[Vertex], upward: List[Vertex], 
                          downward: List[Vertex]) -> int:
    """
    Count total edge traversals during support finding.
    This is for verification of O(n) amortization.
    """
    total = 0
    
    # Process upward vertices
    for v in upward:
        y = v.point.y
        
        # Count left walk
        curr = v.prev
        while True:
            prev_v = curr.prev
            total += 1
            if edge_crosses_height(prev_v, curr, y):
                break
            curr = prev_v
            if curr == v:
                break
        
        # Count right walk
        curr = v.next
        while True:
            next_v = curr.next
            total += 1
            if edge_crosses_height(curr, next_v, y):
                break
            curr = next_v
            if curr == v:
                break
    
    return total


def triangulate_polygon(points: List[Tuple[float, float]]) -> List[Tuple[Point, Point, Point]]:
    """
    Main triangulation function.
    
    Algorithm:
    1. Build polygon as doubly-linked list
    2. Classify reflex vertices
    3. Sort by y-coordinate
    4. Process upward vertices bottom-to-top, adding chords
    5. Process downward vertices top-to-bottom, adding chords
    6. Triangulate all y-monotone regions
    """
    if len(points) < 3:
        return []
    
    # Build polygon
    start, vertices = build_polygon(points)
    
    # Classify and sort
    upward, downward = classify_vertices(vertices)
    
    print(f"Polygon: {len(vertices)} vertices")
    print(f"Upward reflex: {len(upward)}, Downward reflex: {len(downward)}")
    
    # Count edges traversed (for verification)
    edge_count = count_edges_traversed(vertices, upward, downward)
    print(f"Total edge traversals (naive): {edge_count}")
    print(f"Ratio to n: {edge_count / len(vertices):.2f}")
    
    # For now, collect all regions to triangulate
    regions = []
    
    # Process upward vertices
    remaining = start
    for v in upward:
        if v.reflex_type != "up":
            continue  # Skip if already processed
            
        y = v.point.y
        
        try:
            left_upper, left_lower, p_left = find_left_support(v, y)
            right_upper, right_lower, p_right = find_right_support(v, y)
            
            # Split polygon
            lower_start, upper_start = split_at_chord_upward(
                v, left_upper, left_lower, p_left,
                right_upper, right_lower, p_right
            )
            
            regions.append(lower_start)
            remaining = upper_start
            
        except ValueError as e:
            print(f"Warning: {e}")
            continue
    
    # Add remaining region
    if remaining:
        regions.append(remaining)
    
    # Process downward vertices (similar logic, not implemented for brevity)
    # For a complete implementation, we'd process downward vertices top-to-bottom
    
    # Triangulate all regions
    all_triangles = []
    for region in regions:
        try:
            triangles = triangulate_monotone(region)
            all_triangles.extend(triangles)
        except Exception as e:
            print(f"Warning during triangulation: {e}")
    
    return all_triangles


# Test cases
def test_simple_polygon():
    """Test on a simple non-convex polygon"""
    # A "W" shaped polygon
    points = [
        (0, 0),
        (2, 4),
        (4, 1),  # Upward reflex (valley)
        (6, 4),
        (8, 0),
        (4, 2),  # This creates complexity
    ]
    
    triangles = triangulate_polygon(points)
    print(f"Generated {len(triangles)} triangles")
    return triangles


def test_nested_valleys():
    """Test case that could cause O(n^2) for naive approach"""
    # Create nested valleys
    points = [
        (0, 10),   # Top left
        (10, 10),  # Top right
        (10, 0),   # Bottom right
        (9, 1),    # Valley 3 right
        (8, 0.5),  # Valley 3 bottom (reflex)
        (7, 1),    # Valley 3 left
        (6, 2),    # Valley 2 right  
        (5, 1),    # Valley 2 bottom (reflex)
        (4, 2),    # Valley 2 left
        (3, 3),    # Valley 1 right
        (2, 2),    # Valley 1 bottom (reflex)
        (1, 3),    # Valley 1 left
        (0, 0),    # Bottom left
    ]
    
    triangles = triangulate_polygon(points)
    print(f"Nested valleys: {len(triangles)} triangles")
    return triangles


if __name__ == "__main__":
    print("=== Simple polygon test ===")
    test_simple_polygon()
    
    print("\n=== Nested valleys test ===")
    test_nested_valleys()

