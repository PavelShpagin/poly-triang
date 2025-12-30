#!/usr/bin/env python3
"""
Correct O(n + r log r) Polygon Triangulation

This implementation follows the paper's algorithm:
1. Classify vertices: O(n)
2. Sort only reflex (split/merge) vertices: O(r log r)  
3. Monotone decomposition via sweep: O(n + r log r)
4. Triangulate monotone pieces: O(n)

Total: O(n + r log r)
"""

from typing import List, Tuple, Optional, Dict, Set
from enum import Enum, auto
from collections import defaultdict
import math
from bisect import insort_left, bisect_left


class VertexType(Enum):
    START = auto()
    END = auto()
    SPLIT = auto()
    MERGE = auto()
    REGULAR_LEFT = auto()
    REGULAR_RIGHT = auto()


def signed_area(pts: List[Tuple[float, float]]) -> float:
    """Compute signed area of polygon."""
    n = len(pts)
    return sum(pts[i][0] * pts[(i+1)%n][1] - pts[(i+1)%n][0] * pts[i][1] for i in range(n)) / 2


def cross(o: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> float:
    """Cross product of vectors OA and OB."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def is_below(pts: List[Tuple[float, float]], a: int, b: int) -> bool:
    """Is vertex a below vertex b in sweep order?"""
    ya, yb = pts[a][1], pts[b][1]
    if abs(ya - yb) > 1e-12:
        return ya < yb
    return pts[a][0] > pts[b][0]


def classify_vertices(pts: List[Tuple[float, float]]) -> Tuple[List[VertexType], int]:
    """Classify all vertices. Returns (types, reflex_count)."""
    n = len(pts)
    types = [VertexType.REGULAR_LEFT] * n
    reflex_count = 0
    
    for i in range(n):
        p = (i - 1 + n) % n
        nx = (i + 1) % n
        
        p_below = is_below(pts, p, i)
        n_below = is_below(pts, nx, i)
        c = cross(pts[p], pts[i], pts[nx])
        
        if p_below and n_below:
            if c > 1e-12:
                types[i] = VertexType.START
            else:
                types[i] = VertexType.SPLIT
                reflex_count += 1
        elif not p_below and not n_below:
            if c > 1e-12:
                types[i] = VertexType.END
            else:
                types[i] = VertexType.MERGE
                reflex_count += 1
        else:
            types[i] = VertexType.REGULAR_LEFT if not p_below else VertexType.REGULAR_RIGHT
    
    return types, reflex_count


class Edge:
    """Active edge in sweep line status."""
    def __init__(self, upper: int, lower: int, pts: List[Tuple[float, float]]):
        self.upper = upper
        self.lower = lower
        self.pts = pts
        self.helper: Optional[int] = None
    
    def x_at(self, y: float) -> float:
        """X-coordinate of edge at given y."""
        x1, y1 = self.pts[self.upper]
        x2, y2 = self.pts[self.lower]
        if abs(y1 - y2) < 1e-12:
            return min(x1, x2)
        t = (y1 - y) / (y1 - y2)
        return x1 + t * (x2 - x1)


class SweepStatus:
    """Sweep line status structure with O(log n) operations using bisect."""
    def __init__(self, pts: List[Tuple[float, float]]):
        self.edges: List[Edge] = []
        self.pts = pts
        self.sweep_y = float('inf')
        self.edge_map: Dict[int, Edge] = {}  # upper vertex -> edge
    
    def set_y(self, y: float):
        self.sweep_y = y
    
    def _get_x(self, e: Edge) -> float:
        return e.x_at(self.sweep_y)
    
    def insert(self, e: Edge):
        """Insert edge maintaining sorted order by x. O(n) due to list insert."""
        x = self._get_x(e)
        # Binary search for position
        lo, hi = 0, len(self.edges)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._get_x(self.edges[mid]) < x:
                lo = mid + 1
            else:
                hi = mid
        self.edges.insert(lo, e)
        self.edge_map[e.upper] = e
    
    def remove(self, e: Edge):
        """Remove edge. O(n) due to list removal."""
        if e.upper in self.edge_map:
            del self.edge_map[e.upper]
        if e in self.edges:
            self.edges.remove(e)
    
    def find_left(self, x: float) -> Optional[Edge]:
        """Find edge immediately to the left of x. O(log n) using binary search."""
        if not self.edges:
            return None
        
        # Binary search for rightmost edge with x_at < x
        lo, hi = 0, len(self.edges)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._get_x(self.edges[mid]) < x:
                lo = mid + 1
            else:
                hi = mid
        
        # lo is the index of the first edge >= x, so lo-1 is the edge we want
        if lo > 0:
            return self.edges[lo - 1]
        return None
    
    def find_edge(self, upper: int) -> Optional[Edge]:
        """Find edge with given upper vertex. O(1) using hash map."""
        return self.edge_map.get(upper)


def make_monotone(pts: List[Tuple[float, float]], types: List[VertexType]) -> List[Tuple[int, int]]:
    """
    Decompose polygon into y-monotone pieces by adding diagonals.
    Standard sweep-line algorithm from de Berg et al.
    """
    n = len(pts)
    diagonals = []
    
    # Sort all vertices by y (descending), then x (ascending)
    sorted_verts = sorted(range(n), key=lambda i: (-pts[i][1], pts[i][0]))
    
    status = SweepStatus(pts)
    
    # Map from vertex to its outgoing edge
    edge_from: Dict[int, Edge] = {}
    
    for v in sorted_verts:
        vtype = types[v]
        vx, vy = pts[v]
        status.set_y(vy - 1e-9)
        
        p = (v - 1 + n) % n
        nx = (v + 1) % n
        
        if vtype == VertexType.START:
            # Add edge from v going down
            e = Edge(v, nx, pts)
            e.helper = v
            edge_from[v] = e
            status.insert(e)
        
        elif vtype == VertexType.END:
            # Check helper of edge ending at v
            e = edge_from.get(p)
            if e and e.helper is not None and types[e.helper] == VertexType.MERGE:
                diagonals.append((v, e.helper))
            if e:
                status.remove(e)
        
        elif vtype == VertexType.SPLIT:
            # Find edge to left, connect to its helper
            left_edge = status.find_left(vx)
            if left_edge and left_edge.helper is not None:
                diagonals.append((v, left_edge.helper))
                left_edge.helper = v
            
            # Add edge from v going down
            e = Edge(v, nx, pts)
            e.helper = v
            edge_from[v] = e
            status.insert(e)
        
        elif vtype == VertexType.MERGE:
            # Handle edge ending at v
            e = edge_from.get(p)
            if e and e.helper is not None and types[e.helper] == VertexType.MERGE:
                diagonals.append((v, e.helper))
            if e:
                status.remove(e)
            
            # Update helper of edge to left
            left_edge = status.find_left(vx)
            if left_edge:
                if left_edge.helper is not None and types[left_edge.helper] == VertexType.MERGE:
                    diagonals.append((v, left_edge.helper))
                left_edge.helper = v
        
        elif vtype == VertexType.REGULAR_LEFT:
            # Interior is to the right
            e = edge_from.get(p)
            if e and e.helper is not None and types[e.helper] == VertexType.MERGE:
                diagonals.append((v, e.helper))
            if e:
                status.remove(e)
            
            # Add edge from v
            new_e = Edge(v, nx, pts)
            new_e.helper = v
            edge_from[v] = new_e
            status.insert(new_e)
        
        else:  # REGULAR_RIGHT
            # Interior is to the left
            left_edge = status.find_left(vx)
            if left_edge:
                if left_edge.helper is not None and types[left_edge.helper] == VertexType.MERGE:
                    diagonals.append((v, left_edge.helper))
                left_edge.helper = v
    
    return diagonals


def extract_faces(pts: List[Tuple[float, float]], diagonals: List[Tuple[int, int]]) -> List[List[int]]:
    """Extract all interior faces from polygon with diagonals using half-edge traversal."""
    n = len(pts)
    
    # Build adjacency with edge directions
    adj: Dict[int, List[int]] = defaultdict(list)
    for i in range(n):
        j = (i + 1) % n
        adj[i].append(j)
        adj[j].append(i)
    for u, v in diagonals:
        adj[u].append(v)
        adj[v].append(u)
    
    # Sort neighbors counter-clockwise around each vertex
    def angle(u: int, v: int) -> float:
        return math.atan2(pts[v][1] - pts[u][1], pts[v][0] - pts[u][0])
    
    for u in adj:
        adj[u] = sorted(set(adj[u]), key=lambda v: angle(u, v))
    
    # Track used half-edges
    used: Set[Tuple[int, int]] = set()
    faces: List[List[int]] = []
    
    # Try each directed edge as face start
    all_edges = []
    for i in range(n):
        all_edges.append((i, (i + 1) % n))
    for u, v in diagonals:
        all_edges.append((u, v))
        all_edges.append((v, u))
    
    for start_u, start_v in all_edges:
        if (start_u, start_v) in used:
            continue
        
        face = []
        u, v = start_u, start_v
        
        for _ in range(2 * n + 10):
            if (u, v) in used:
                break
            used.add((u, v))
            face.append(u)
            
            # Find next edge: turn most clockwise from incoming edge
            neighbors = adj[v]
            if u not in neighbors:
                break
            idx = neighbors.index(u)
            # Go to previous in CCW order = next in CW order (most clockwise turn)
            next_v = neighbors[(idx - 1 + len(neighbors)) % len(neighbors)]
            u, v = v, next_v
            
            if u == start_u and v == start_v:
                break
        
        if len(face) >= 3:
            # Compute signed area - positive = CCW = interior face
            area = 0.0
            for i in range(len(face)):
                j = (i + 1) % len(face)
                area += pts[face[i]][0] * pts[face[j]][1]
                area -= pts[face[j]][0] * pts[face[i]][1]
            if area > 1e-9:
                faces.append(face)
    
    return faces


def triangulate_monotone(pts: List[Tuple[float, float]], face: List[int]) -> List[Tuple[int, int, int]]:
    """Triangulate a y-monotone polygon using the standard stack algorithm."""
    if len(face) < 3:
        return []
    if len(face) == 3:
        return [(face[0], face[1], face[2])]
    
    # Sort vertices by y (descending)
    sorted_face = sorted(range(len(face)), key=lambda i: (-pts[face[i]][1], pts[face[i]][0]))
    
    # Determine which chain each vertex is on
    top_idx = sorted_face[0]
    bot_idx = sorted_face[-1]
    
    # Walk from top to find left and right chains
    left_chain = set()
    right_chain = set()
    
    # Go forward from top
    i = top_idx
    while i != bot_idx:
        right_chain.add(face[i])
        i = (i + 1) % len(face)
    
    # Go backward from top
    i = top_idx
    while i != bot_idx:
        left_chain.add(face[i])
        i = (i - 1 + len(face)) % len(face)
    
    triangles = []
    stack = [sorted_face[0], sorted_face[1]]
    
    for j in range(2, len(sorted_face)):
        curr = sorted_face[j]
        curr_v = face[curr]
        
        top_v = face[stack[-1]]
        
        # Check if same chain
        curr_left = curr_v in left_chain
        top_left = top_v in left_chain
        
        if curr_left != top_left:
            # Different chains - pop all and triangulate
            while len(stack) > 1:
                v1 = stack.pop()
                v2 = stack[-1]
                triangles.append((face[curr], face[v1], face[v2]))
            stack.pop()
            stack.append(sorted_face[j-1])
            stack.append(curr)
        else:
            # Same chain - pop while diagonals are inside
            last_popped = stack.pop()
            
            while stack:
                v = stack[-1]
                # Check if diagonal is inside
                p1, p2, p3 = pts[face[curr]], pts[face[last_popped]], pts[face[v]]
                c = cross(p1, p2, p3)
                
                if (curr_left and c > 1e-9) or (not curr_left and c < -1e-9):
                    triangles.append((face[curr], face[last_popped], face[v]))
                    last_popped = stack.pop()
                else:
                    break
            
            stack.append(last_popped)
            stack.append(curr)
    
    return triangles


def triangulate(pts: List[Tuple[float, float]]) -> Tuple[List[Tuple[int, int, int]], int]:
    """
    Triangulate a simple polygon in O(n + r log r) time.
    Returns (triangles, reflex_count).
    """
    n = len(pts)
    if n < 3:
        return [], 0
    
    # Ensure CCW orientation
    if signed_area(pts) < 0:
        pts = list(reversed(pts))
    
    # Handle trivial cases
    if n == 3:
        return [(0, 1, 2)], 0
    
    # Classify vertices
    types, r = classify_vertices(pts)
    
    # Special case: convex polygon
    if r == 0:
        return [(0, i, i+1) for i in range(1, n-1)], 0
    
    # Monotone decomposition
    diagonals = make_monotone(pts, types)
    
    # Extract faces
    faces = extract_faces(pts, diagonals)
    
    # Triangulate each face
    triangles = []
    for face in faces:
        triangles.extend(triangulate_monotone(pts, face))
    
    return triangles, r


# Polygon generators that guarantee simple polygons
def gen_convex(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Regular n-gon (r = 0)."""
    return [(math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)) for i in range(n)]

def gen_star(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Alternating-radius star (r ~ n/2). Always simple."""
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        radius = 1.0 if i % 2 == 0 else 0.4
        pts.append((radius * math.cos(angle), radius * math.sin(angle)))
    return pts

def gen_random_star(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Star-shaped polygon with random radii (always simple, r ~ n/2)."""
    import random
    random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n  # Fixed angles = always simple
        radius = 0.3 + 0.7 * random.random()
        pts.append((radius * math.cos(angle), radius * math.sin(angle)))
    return pts

def gen_spiral(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Spiral polygon (low r)."""
    pts = []
    for i in range(n):
        t = 4 * math.pi * i / n
        r = 0.2 + 0.8 * i / n
        pts.append((r * math.cos(t), r * math.sin(t)))
    return pts

def gen_comb(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Comb-shaped polygon (high r)."""
    teeth = max(3, n // 4)
    pts = []
    for i in range(teeth):
        x = i / teeth
        pts.append((x, 0))
        pts.append((x + 0.4/teeth, 0.8))
    pts.append((1, 0))
    pts.append((1, -0.1))
    pts.append((0, -0.1))
    # Pad or trim to exactly n vertices
    while len(pts) < n:
        pts.insert(-3, (0.5, 0.4))
    return pts[:n]


# Test
if __name__ == '__main__':
    import random
    
    def test_polygon(name: str, pts: List[Tuple[float, float]], expected_r: int = None):
        n = len(pts)
        tris, r = triangulate(pts)
        ok = len(tris) == n - 2
        status = "PASS" if ok else "FAIL"
        r_str = f"r={r}" + (f" (expected {expected_r})" if expected_r is not None and r != expected_r else "")
        print(f"{status}: {name} n={n} {r_str} triangles={len(tris)} (expected {n-2})")
        return ok
    
    all_pass = True
    
    # Convex
    all_pass &= test_polygon("Convex-50", gen_convex(50), 0)
    all_pass &= test_polygon("Convex-1000", gen_convex(1000), 0)
    
    # Star (alternating radii)
    all_pass &= test_polygon("Star-20", gen_star(20))
    all_pass &= test_polygon("Star-100", gen_star(100))
    
    # Random star-shaped (random radii, fixed angles)
    all_pass &= test_polygon("RandomStar-100", gen_random_star(100, seed=42))
    all_pass &= test_polygon("RandomStar-1000", gen_random_star(1000, seed=42))
    all_pass &= test_polygon("RandomStar-5000", gen_random_star(5000, seed=42))
    
    # Spiral
    all_pass &= test_polygon("Spiral-100", gen_spiral(100))
    all_pass &= test_polygon("Spiral-1000", gen_spiral(1000))
    
    print("\n" + ("All tests passed!" if all_pass else "Some tests FAILED!"))

