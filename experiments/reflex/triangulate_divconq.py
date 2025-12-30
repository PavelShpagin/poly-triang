"""
Divide and Conquer Triangulation

Key idea: Find a diagonal that roughly balances reflex vertices,
then recurse.

If we can find such a diagonal in O(k) time where k = vertices in subproblem,
we get T(n,r) = T(n1,r/2) + T(n2,r/2) + O(n) = O(n log r)

Better: If diagonal finding is O(r), we get O(n + r log r).
"""

import math
from typing import List, Tuple, Optional, Set
from dataclasses import dataclass


@dataclass
class Vertex:
    x: float
    y: float
    idx: int
    is_reflex: bool = False


def cross2d(ox, oy, ax, ay, bx, by):
    return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox)


def ccw(a, b, c):
    return cross2d(a.x, a.y, b.x, b.y, c.x, c.y) > 0


def segments_intersect(a1, a2, b1, b2):
    """Check if segment a1-a2 intersects segment b1-b2 (not at endpoints)."""
    d1 = cross2d(b1.x, b1.y, b2.x, b2.y, a1.x, a1.y)
    d2 = cross2d(b1.x, b1.y, b2.x, b2.y, a2.x, a2.y)
    d3 = cross2d(a1.x, a1.y, a2.x, a2.y, b1.x, b1.y)
    d4 = cross2d(a1.x, a1.y, a2.x, a2.y, b2.x, b2.y)
    
    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and \
       ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    return False


def point_in_triangle(p, a, b, c):
    """Check if p is strictly inside triangle abc."""
    d1 = cross2d(a.x, a.y, b.x, b.y, p.x, p.y)
    d2 = cross2d(b.x, b.y, c.x, c.y, p.x, p.y)
    d3 = cross2d(c.x, c.y, a.x, a.y, p.x, p.y)
    
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    
    return not (has_neg and has_pos)


class DivideConquerTriangulator:
    """
    Divide and conquer triangulation.
    
    1. Find a diagonal that splits reflex vertices roughly in half
    2. Recurse on both subpolygons
    3. Merge results
    """
    
    def __init__(self, points: List[Tuple[float, float]]):
        self.original_points = points
        self.triangles = []
        
        self.stats = {
            'n': len(points),
            'r': 0,
            'diagonal_searches': 0,
            'recursion_depth': 0,
            'max_depth': 0,
        }
    
    def is_convex_vertex(self, verts: List[Vertex], i: int) -> bool:
        """Check if vertex i is convex in the polygon."""
        n = len(verts)
        prev_v = verts[(i - 1) % n]
        curr_v = verts[i]
        next_v = verts[(i + 1) % n]
        return cross2d(prev_v.x, prev_v.y, curr_v.x, curr_v.y, next_v.x, next_v.y) > 0
    
    def diagonal_valid(self, verts: List[Vertex], i: int, j: int) -> bool:
        """Check if diagonal from vertex i to vertex j is valid (inside polygon, no intersections)."""
        n = len(verts)
        if abs(i - j) <= 1 or abs(i - j) == n - 1:
            return False  # Adjacent vertices
        
        vi, vj = verts[i], verts[j]
        
        # Check if diagonal is inside the polygon at both endpoints
        # At vertex i: diagonal should be on the interior side
        prev_i = verts[(i - 1) % n]
        next_i = verts[(i + 1) % n]
        
        # For a valid diagonal from i, vj should be in the "cone" from prev_i to next_i
        # This is complex, let's use a simpler check
        
        # Check if midpoint of diagonal is inside the polygon (approximate)
        mid_x = (vi.x + vj.x) / 2
        mid_y = (vi.y + vj.y) / 2
        
        # Check intersection with all edges
        for k in range(n):
            k_next = (k + 1) % n
            if k == i or k_next == i or k == j or k_next == j:
                continue
            
            if segments_intersect(vi, vj, verts[k], verts[k_next]):
                return False
        
        return True
    
    def find_splitting_diagonal(self, verts: List[Vertex], reflex_indices: Set[int]) -> Optional[Tuple[int, int]]:
        """
        Find a diagonal that splits reflex vertices roughly in half.
        
        Strategy: Try diagonals from reflex vertices first.
        """
        n = len(verts)
        r = len(reflex_indices)
        
        if r <= 1:
            # Any valid diagonal works
            for i in range(n):
                for j in range(i + 2, n):
                    if j == (i + n - 1) % n:
                        continue
                    if self.diagonal_valid(verts, i, j):
                        self.stats['diagonal_searches'] += 1
                        return (i, j)
            return None
        
        # Try diagonals from each reflex vertex
        reflex_list = list(reflex_indices)
        target = r // 2
        
        best_diagonal = None
        best_balance = float('inf')
        
        for ri in reflex_list[:min(10, len(reflex_list))]:  # Limit search
            for j in range(n):
                if abs(ri - j) <= 1 or abs(ri - j) == n - 1:
                    continue
                
                self.stats['diagonal_searches'] += 1
                
                if not self.diagonal_valid(verts, ri, j):
                    continue
                
                # Count reflex vertices on each side
                # Side 1: from ri to j (exclusive)
                # Side 2: from j to ri (exclusive)
                
                side1_reflex = 0
                side2_reflex = 0
                
                k = (ri + 1) % n
                while k != j:
                    if k in reflex_indices:
                        side1_reflex += 1
                    k = (k + 1) % n
                
                k = (j + 1) % n
                while k != ri:
                    if k in reflex_indices:
                        side2_reflex += 1
                    k = (k + 1) % n
                
                balance = abs(side1_reflex - side2_reflex)
                if balance < best_balance:
                    best_balance = balance
                    best_diagonal = (ri, j)
                
                if best_balance == 0:
                    return best_diagonal
        
        return best_diagonal
    
    def triangulate_convex(self, verts: List[Vertex]) -> List[Tuple[int, int, int]]:
        """Fan triangulate a convex polygon."""
        tris = []
        n = len(verts)
        for i in range(1, n - 1):
            tris.append((verts[0].idx, verts[i].idx, verts[i + 1].idx))
        return tris
    
    def triangulate_recursive(self, verts: List[Vertex], depth: int) -> List[Tuple[int, int, int]]:
        """Recursive triangulation."""
        self.stats['recursion_depth'] = depth
        self.stats['max_depth'] = max(self.stats['max_depth'], depth)
        
        n = len(verts)
        
        if n < 3:
            return []
        
        if n == 3:
            return [(verts[0].idx, verts[1].idx, verts[2].idx)]
        
        # Classify vertices
        reflex_indices = set()
        for i in range(n):
            if not self.is_convex_vertex(verts, i):
                reflex_indices.add(i)
        
        r = len(reflex_indices)
        
        if r == 0:
            # Convex polygon - fan triangulate
            return self.triangulate_convex(verts)
        
        # Find a splitting diagonal
        diagonal = self.find_splitting_diagonal(verts, reflex_indices)
        
        if diagonal is None:
            # Fallback: ear clipping
            return self.ear_clip(verts)
        
        i, j = diagonal
        
        # Split into two subpolygons
        # Subpoly 1: i to j (including both)
        # Subpoly 2: j to i (including both)
        
        sub1 = []
        k = i
        while True:
            sub1.append(verts[k])
            if k == j:
                break
            k = (k + 1) % n
        
        sub2 = []
        k = j
        while True:
            sub2.append(verts[k])
            if k == i:
                break
            k = (k + 1) % n
        
        # Recurse
        tris1 = self.triangulate_recursive(sub1, depth + 1)
        tris2 = self.triangulate_recursive(sub2, depth + 1)
        
        return tris1 + tris2
    
    def ear_clip(self, verts: List[Vertex]) -> List[Tuple[int, int, int]]:
        """Simple ear clipping fallback."""
        tris = []
        verts = list(verts)  # Copy
        
        while len(verts) > 3:
            n = len(verts)
            found = False
            
            for i in range(n):
                prev_v = verts[(i - 1) % n]
                curr_v = verts[i]
                next_v = verts[(i + 1) % n]
                
                # Check if convex
                if cross2d(prev_v.x, prev_v.y, curr_v.x, curr_v.y, next_v.x, next_v.y) <= 0:
                    continue
                
                # Check if ear (no other vertex inside)
                is_ear = True
                for j in range(n):
                    if j == (i - 1) % n or j == i or j == (i + 1) % n:
                        continue
                    if point_in_triangle(verts[j], prev_v, curr_v, next_v):
                        is_ear = False
                        break
                
                if is_ear:
                    tris.append((prev_v.idx, curr_v.idx, next_v.idx))
                    verts.pop(i)
                    found = True
                    break
            
            if not found:
                break
        
        if len(verts) == 3:
            tris.append((verts[0].idx, verts[1].idx, verts[2].idx))
        
        return tris
    
    def triangulate(self) -> Tuple[List[Tuple[int, int, int]], dict]:
        verts = [Vertex(x, y, i) for i, (x, y) in enumerate(self.original_points)]
        
        # Count reflex vertices
        n = len(verts)
        for i in range(n):
            if not self.is_convex_vertex(verts, i):
                verts[i].is_reflex = True
                self.stats['r'] += 1
        
        self.triangles = self.triangulate_recursive(verts, 0)
        
        return self.triangles, self.stats


# ==============================================================
# Test
# ==============================================================

def create_star(n):
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n - math.pi/2
        r = 2 if i % 2 == 0 else 1
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def create_convex(n):
    return [(math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)) for i in range(n)]


def create_comb(n):
    pts = []
    teeth = n // 4
    for i in range(teeth):
        pts.append((4*i, 0))
        pts.append((4*i + 1, 2))
        pts.append((4*i + 2, 0))
    pts.append((4*teeth, 0))
    return pts


def benchmark():
    print("=" * 80)
    print("DIVIDE AND CONQUER TRIANGULATION")
    print("=" * 80)
    print()
    
    for name, gen in [("Convex", create_convex), ("Comb", create_comb), ("Star", create_star)]:
        print(f"{name}:")
        print("-" * 70)
        print(f"{'n':>6} {'r':>5} {'searches':>10} {'depth':>8} {'srch/n':>10}")
        print("-" * 70)
        
        for n in [50, 100, 200, 500]:
            pts = gen(n)
            tri = DivideConquerTriangulator(pts)
            _, stats = tri.triangulate()
            
            searches = stats['diagonal_searches']
            depth = stats['max_depth']
            
            print(f"{stats['n']:>6} {stats['r']:>5} {searches:>10} {depth:>8} "
                  f"{searches/stats['n']:>10.2f}")
        print()
    
    print("=" * 80)
    print("ANALYSIS:")
    print("  depth should be O(log r)")
    print("  searches/n ideally O(1) or O(log r)")
    print("=" * 80)


if __name__ == "__main__":
    benchmark()

