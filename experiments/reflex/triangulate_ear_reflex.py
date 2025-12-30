"""
Ear Clipping with Reflex Vertex Tracking

Key insight: A vertex v is a valid ear iff:
1. v is convex (interior angle < 180)
2. No REFLEX vertex is inside triangle(v.prev, v, v.next)

By tracking only reflex vertices, we might achieve O(n + r * something).

The "something" depends on how many triangles each reflex vertex affects.
"""

import math
from typing import List, Tuple, Set, Dict
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class Vertex:
    x: float
    y: float
    idx: int
    prev: 'Vertex' = None
    next: 'Vertex' = None
    is_ear: bool = False
    is_reflex: bool = False


def sign(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


def cross2d(o, a, b):
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def point_in_triangle(p, a, b, c):
    """Check if point p is strictly inside triangle abc."""
    s1 = sign(cross2d(a, b, p))
    s2 = sign(cross2d(b, c, p))
    s3 = sign(cross2d(c, a, p))
    
    # All same sign (and non-zero) means strictly inside
    if s1 == s2 == s3 and s1 != 0:
        return True
    return False


class EarClipperReflexTracking:
    """
    O(n + r * k) ear clipping where k = average updates per reflex vertex.
    
    Key data structures:
    - reflex_vertices: set of reflex vertex indices
    - blocking: for each potential ear, which reflex vertices block it
    - blocked_by: for each reflex vertex, which ears it blocks
    """
    
    def __init__(self, points: List[Tuple[float, float]]):
        self.n = len(points)
        self.vertices = []
        
        # Build circular linked list
        for i, (x, y) in enumerate(points):
            self.vertices.append(Vertex(x, y, i))
        
        for i in range(self.n):
            self.vertices[i].prev = self.vertices[(i - 1) % self.n]
            self.vertices[i].next = self.vertices[(i + 1) % self.n]
        
        self.triangles = []
        self.reflex_set = set()
        
        # Tracking structures
        self.blocking = defaultdict(set)    # vertex_idx -> set of reflex_idx blocking it
        self.blocked_by = defaultdict(set)  # reflex_idx -> set of vertex_idx it blocks
        
        self.stats = {
            'n': self.n,
            'r': 0,
            'ear_checks': 0,
            'blocking_updates': 0,
            'triangles_created': 0,
        }
    
    def is_convex(self, v: Vertex) -> bool:
        """Check if vertex v is convex (interior angle < 180)."""
        return cross2d(v.prev, v, v.next) > 0
    
    def classify_vertices(self):
        """Classify all vertices as convex or reflex: O(n)."""
        for v in self.vertices:
            v.is_reflex = not self.is_convex(v)
            if v.is_reflex:
                self.reflex_set.add(v.idx)
        self.stats['r'] = len(self.reflex_set)
    
    def compute_initial_blocking(self):
        """
        For each convex vertex, compute which reflex vertices block it.
        
        Naive: O(n * r)
        """
        for v in self.vertices:
            if v.is_reflex:
                continue
            
            # Check which reflex vertices are inside triangle(v.prev, v, v.next)
            for ridx in self.reflex_set:
                r = self.vertices[ridx]
                if r == v or r == v.prev or r == v.next:
                    continue
                
                if point_in_triangle(r, v.prev, v, v.next):
                    self.blocking[v.idx].add(ridx)
                    self.blocked_by[ridx].add(v.idx)
                    self.stats['ear_checks'] += 1
            
            # v is an ear if no reflex vertices block it
            v.is_ear = len(self.blocking[v.idx]) == 0
    
    def update_vertex(self, v: Vertex):
        """Update ear status of vertex v after a neighbor was removed."""
        if v.is_reflex:
            # Check if still reflex
            v.is_reflex = not self.is_convex(v)
            if not v.is_reflex:
                # Was reflex, now convex
                self.reflex_set.discard(v.idx)
                # Remove from blocked_by for all ears it was blocking
                for ear_idx in list(self.blocked_by[v.idx]):
                    self.blocking[ear_idx].discard(v.idx)
                    # Check if that ear is now valid
                    if len(self.blocking[ear_idx]) == 0:
                        self.vertices[ear_idx].is_ear = True
                self.blocked_by[v.idx].clear()
                self.stats['blocking_updates'] += 1
        
        # Recompute blocking for v (as potential ear)
        old_blocking = self.blocking[v.idx].copy()
        self.blocking[v.idx].clear()
        
        for ridx in self.reflex_set:
            r = self.vertices[ridx]
            if r == v or r == v.prev or r == v.next:
                continue
            
            if point_in_triangle(r, v.prev, v, v.next):
                self.blocking[v.idx].add(ridx)
                self.blocked_by[ridx].add(v.idx)
            else:
                self.blocked_by[ridx].discard(v.idx)
            
            self.stats['ear_checks'] += 1
        
        v.is_ear = not v.is_reflex and len(self.blocking[v.idx]) == 0
        self.stats['blocking_updates'] += 1
    
    def clip_ear(self, v: Vertex):
        """Remove ear v and form a triangle."""
        # Create triangle
        self.triangles.append((v.prev.idx, v.idx, v.next.idx))
        self.stats['triangles_created'] += 1
        
        # Remove v from linked list
        v.prev.next = v.next
        v.next.prev = v.prev
        
        # Clear blocking info for v
        for ridx in self.blocking[v.idx]:
            self.blocked_by[ridx].discard(v.idx)
        self.blocking[v.idx].clear()
        
        # If v was reflex (shouldn't happen, but just in case)
        if v.idx in self.reflex_set:
            self.reflex_set.discard(v.idx)
            for ear_idx in self.blocked_by[v.idx]:
                self.blocking[ear_idx].discard(v.idx)
        
        # Update neighbors
        self.update_vertex(v.prev)
        self.update_vertex(v.next)
    
    def triangulate(self) -> Tuple[List[Tuple[int, int, int]], Dict]:
        """Main triangulation routine."""
        self.classify_vertices()
        self.compute_initial_blocking()
        
        remaining = self.n
        current = self.vertices[0]
        
        # Keep clipping ears until only 3 vertices remain
        iterations = 0
        max_iterations = self.n * self.n  # Safety limit
        
        while remaining > 3 and iterations < max_iterations:
            iterations += 1
            
            # Find an ear
            found_ear = False
            start = current
            
            while True:
                if current.is_ear:
                    found_ear = True
                    break
                current = current.next
                if current == start:
                    break
            
            if not found_ear:
                # No ear found - this shouldn't happen for a simple polygon
                # Try to force progress by picking any convex vertex
                current = start
                while True:
                    if not current.is_reflex:
                        found_ear = True
                        break
                    current = current.next
                    if current == start:
                        break
            
            if found_ear:
                next_v = current.next
                self.clip_ear(current)
                remaining -= 1
                current = next_v
            else:
                # Stuck - polygon might not be simple
                break
        
        # Final triangle
        if remaining == 3:
            self.triangles.append((current.prev.idx, current.idx, current.next.idx))
            self.stats['triangles_created'] += 1
        
        return self.triangles, self.stats


# ==============================================================
# Alternative: Smarter Blocking Updates
# ==============================================================

class SmartEarClipper:
    """
    Improved ear clipping with smarter blocking updates.
    
    Key optimization: When an ear is clipped, only check reflex vertices
    that were "near" the clipped triangle.
    """
    
    def __init__(self, points: List[Tuple[float, float]]):
        self.n = len(points)
        self.vertices = []
        
        for i, (x, y) in enumerate(points):
            self.vertices.append(Vertex(x, y, i))
        
        for i in range(self.n):
            self.vertices[i].prev = self.vertices[(i - 1) % self.n]
            self.vertices[i].next = self.vertices[(i + 1) % self.n]
        
        self.triangles = []
        self.reflex_set = set()
        self.ear_list = []  # List of ear vertices
        
        self.stats = {
            'n': self.n,
            'r': 0,
            'point_in_tri_tests': 0,
            'ear_updates': 0,
        }
    
    def is_convex(self, v: Vertex) -> bool:
        return cross2d(v.prev, v, v.next) > 0
    
    def is_ear(self, v: Vertex) -> bool:
        """Check if v is a valid ear (convex + no reflex inside)."""
        if v.is_reflex:
            return False
        
        # Check only reflex vertices
        for ridx in self.reflex_set:
            r = self.vertices[ridx]
            if r == v or r == v.prev or r == v.next:
                continue
            
            self.stats['point_in_tri_tests'] += 1
            if point_in_triangle(r, v.prev, v, v.next):
                return False
        
        return True
    
    def triangulate(self) -> Tuple[List[Tuple[int, int, int]], Dict]:
        # Classify vertices
        for v in self.vertices:
            v.is_reflex = not self.is_convex(v)
            if v.is_reflex:
                self.reflex_set.add(v.idx)
        
        self.stats['r'] = len(self.reflex_set)
        
        # Find initial ears
        for v in self.vertices:
            if self.is_ear(v):
                v.is_ear = True
                self.ear_list.append(v)
        
        remaining = self.n
        
        while remaining > 3 and self.ear_list:
            # Get an ear
            v = self.ear_list.pop(0)
            
            # Check if still valid (might have changed)
            if v.prev.next != v or not v.is_ear:
                continue
            
            # Clip the ear
            self.triangles.append((v.prev.idx, v.idx, v.next.idx))
            
            # Update linked list
            prev_v, next_v = v.prev, v.next
            prev_v.next = next_v
            next_v.prev = prev_v
            
            # Remove v from reflex set if it was there
            self.reflex_set.discard(v.idx)
            
            remaining -= 1
            
            # Update prev_v
            prev_v.is_reflex = not self.is_convex(prev_v)
            if not prev_v.is_reflex:
                self.reflex_set.discard(prev_v.idx)
                if self.is_ear(prev_v):
                    if not prev_v.is_ear:
                        prev_v.is_ear = True
                        self.ear_list.append(prev_v)
                else:
                    prev_v.is_ear = False
            self.stats['ear_updates'] += 1
            
            # Update next_v
            next_v.is_reflex = not self.is_convex(next_v)
            if not next_v.is_reflex:
                self.reflex_set.discard(next_v.idx)
                if self.is_ear(next_v):
                    if not next_v.is_ear:
                        next_v.is_ear = True
                        self.ear_list.append(next_v)
                else:
                    next_v.is_ear = False
            self.stats['ear_updates'] += 1
        
        # Final triangle
        if remaining == 3:
            v = self.vertices[0]
            while v.prev.next != v:
                v = self.vertices[v.idx + 1] if v.idx + 1 < self.n else self.vertices[0]
            self.triangles.append((v.prev.idx, v.idx, v.next.idx))
        
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
    """Comb polygon with many teeth."""
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
    print("EAR CLIPPING WITH REFLEX TRACKING")
    print("=" * 80)
    print()
    
    for name, gen in [("Convex", create_convex), ("Comb", create_comb), ("Star", create_star)]:
        print(f"{name}:")
        print("-" * 75)
        print(f"{'n':>6} {'r':>5} {'tests':>10} {'updates':>10} {'tests/n':>10} {'tests/r':>10}")
        print("-" * 75)
        
        for n in [50, 100, 200, 500, 1000]:
            pts = gen(n)
            tri = SmartEarClipper(pts)
            _, stats = tri.triangulate()
            
            r = max(stats['r'], 1)
            tests = stats['point_in_tri_tests']
            updates = stats['ear_updates']
            
            print(f"{stats['n']:>6} {stats['r']:>5} {tests:>10} {updates:>10} "
                  f"{tests/stats['n']:>10.1f} {tests/r:>10.1f}")
        print()
    
    print("=" * 80)
    print("ANALYSIS:")
    print("  If tests/n = O(1), algorithm is O(n)")
    print("  If tests/r = O(n), algorithm is O(nr) = O(n^2) for star")
    print("=" * 80)


if __name__ == "__main__":
    benchmark()

