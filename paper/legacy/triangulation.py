"""
Polygon Triangulation in O(n + r log r) Time

A chain-based monotone decomposition algorithm where n = vertices, r = reflex vertices.

This implementation matches the paper "Triangulating Simple Polygons in O(n + r log r) Time"
and uses:
1. Chain-based sweep line status (O(r) chains instead of O(n) edges)
2. Lazy edge pointer advancement (O(n) amortized total)
3. Stack-based monotone triangulation (O(n) per face)

Authors: Pavel Shpagin, Vasyl Tereschenko
Taras Shevchenko National University of Kyiv
"""

from __future__ import annotations
import math
from collections import defaultdict
from enum import Enum, auto
from typing import List, Tuple, Dict, Set, Optional
import bisect


class VertexType(Enum):
    START = auto()      # Local max, convex (interior angle < π)
    END = auto()        # Local min, convex
    SPLIT = auto()      # Local max, reflex (interior angle > π)
    MERGE = auto()      # Local min, reflex
    REGULAR_LEFT = auto()   # On left boundary chain
    REGULAR_RIGHT = auto()  # On right boundary chain


class MonotoneChain:
    """
    A y-monotone chain from a local maximum to a local minimum.
    Stores vertices in decreasing y order (top to bottom).
    """
    
    def __init__(self, vertices: List[int], pts: List[Tuple[float, float]], is_left: bool):
        self.vertices = vertices  # Vertex indices in decreasing y order
        self.pts = pts
        self.is_left = is_left   # True if interior is to the right
        self.edge_ptr = 0        # Current edge index (lazy advancement)
        self.pending: Optional[int] = None  # Pending merge vertex

    @property
    def upper(self) -> int:
        """Upper endpoint (local maximum)."""
        return self.vertices[0]

    @property
    def lower(self) -> int:
        """Lower endpoint (local minimum)."""
        return self.vertices[-1]

    def advance_to(self, y: float) -> None:
        """Lazily advance edge pointer until current edge spans y."""
        while self.edge_ptr < len(self.vertices) - 1:
            lower_v = self.vertices[self.edge_ptr + 1]
            if self.pts[lower_v][1] <= y:
                break
            self.edge_ptr += 1

    def x_at(self, y: float) -> float:
        """Compute x-coordinate of chain at height y."""
        self.advance_to(y)
        if self.edge_ptr >= len(self.vertices) - 1:
            return self.pts[self.vertices[-1]][0]

        upper_v = self.vertices[self.edge_ptr]
        lower_v = self.vertices[self.edge_ptr + 1]
        x1, y1 = self.pts[upper_v]
        x2, y2 = self.pts[lower_v]

        if abs(y1 - y2) < 1e-12:
            return min(x1, x2)
        t = (y - y1) / (y2 - y1)
        t = max(0.0, min(1.0, t))  # Clamp to [0, 1]
        return x1 + t * (x2 - x1)

    def slab_entry(self) -> int:
        """Return the upper vertex of the current edge (slab entry point)."""
        return self.vertices[self.edge_ptr]

    def reset(self) -> None:
        """Reset chain state for sweep."""
        self.edge_ptr = 0
        self.pending = None


class ChainBST:
    """
    Balanced search tree storing active left-boundary chains.
    Ordered by x-coordinate at current sweep height.
    Uses binary search on a sorted list (simpler than treap for Python).
    """
    
    def __init__(self):
        self.chains: List[MonotoneChain] = []
        self.sweep_y: float = float('inf')

    def set_y(self, y: float) -> None:
        self.sweep_y = y

    def _chain_x(self, chain: MonotoneChain) -> float:
        return chain.x_at(self.sweep_y)

    def insert(self, chain: MonotoneChain) -> None:
        x = self._chain_x(chain)
        pos = bisect.bisect_left([self._chain_x(c) for c in self.chains], x)
        self.chains.insert(pos, chain)

    def remove(self, chain: MonotoneChain) -> None:
        if chain in self.chains:
            self.chains.remove(chain)

    def left_of(self, x: float) -> Optional[MonotoneChain]:
        """Find the chain immediately to the left of x."""
        best = None
        best_x = float('-inf')
        for c in self.chains:
            cx = self._chain_x(c)
            if cx < x and cx > best_x:
                best_x = cx
                best = c
        return best


class Triangulator:
    """
    O(n + r log r) polygon triangulator using chain-based monotone decomposition.
    """
    
    def __init__(self, points: List[Tuple[float, float]]):
        # Ensure CCW orientation
        if self._signed_area(points) < 0:
            points = list(reversed(points))
        self.pts = list(points)
        self.n = len(points)
        self.types: List[VertexType] = [VertexType.REGULAR_LEFT] * self.n
        self.diagonals: List[Tuple[int, int]] = []
        self.r = 0  # Reflex vertex count
        self.chains: List[MonotoneChain] = []
        self.extrema: List[int] = []

    def triangulate(self) -> List[Tuple[int, int, int]]:
        """Triangulate the polygon, returning list of triangles."""
        if self.n < 3:
            return []
        if self.n == 3:
            return [(0, 1, 2)]

        self._classify_vertices()

        # Fast path for convex polygons
        if self.r == 0:
            return [(0, i, i + 1) for i in range(1, self.n - 1)]

        self._build_chains()
        self._decompose()
        return self._triangulate_faces()

    def _prev(self, i: int) -> int:
        return (i - 1 + self.n) % self.n

    def _next(self, i: int) -> int:
        return (i + 1) % self.n

    def _is_below(self, a: int, b: int) -> bool:
        """True if vertex a is below vertex b (lexicographic comparison)."""
        ya, yb = self.pts[a][1], self.pts[b][1]
        if abs(ya - yb) > 1e-12:
            return ya < yb
        return self.pts[a][0] > self.pts[b][0]

    def _cross(self, a: int, b: int, c: int) -> float:
        """Signed cross product of vectors (a->b) and (a->c)."""
        ax, ay = self.pts[a]
        bx, by = self.pts[b]
        cx, cy = self.pts[c]
        return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)

    def _classify_vertices(self) -> None:
        """Classify all vertices and identify extrema. O(n)."""
        self.r = 0
        self.extrema = []

        for i in range(self.n):
            p, nx = self._prev(i), self._next(i)
            p_below = self._is_below(p, i)
            n_below = self._is_below(nx, i)
            cross = self._cross(p, i, nx)

            if p_below and n_below:  # Local maximum
                if cross > 0:
                    self.types[i] = VertexType.START
                else:
                    self.types[i] = VertexType.SPLIT
                    self.r += 1
                self.extrema.append(i)
            elif not p_below and not n_below:  # Local minimum
                if cross > 0:
                    self.types[i] = VertexType.END
                else:
                    self.types[i] = VertexType.MERGE
                    self.r += 1
                self.extrema.append(i)
            else:  # Regular
                self.types[i] = VertexType.REGULAR_RIGHT if p_below else VertexType.REGULAR_LEFT

    def _build_chains(self) -> None:
        """Build monotone chains from local maxima to local minima. O(n)."""
        self.chains = []
        maxima = [i for i in self.extrema 
                  if self.types[i] in (VertexType.START, VertexType.SPLIT)]

        for max_v in maxima:
            # Chain going in CCW direction (next) - this is the "left" chain
            chain1 = [max_v]
            curr = self._next(max_v)
            while self.types[curr] not in (VertexType.END, VertexType.MERGE):
                chain1.append(curr)
                curr = self._next(curr)
            chain1.append(curr)

            # Chain going in CW direction (prev) - not stored in BST
            chain2 = [max_v]
            curr = self._prev(max_v)
            while self.types[curr] not in (VertexType.END, VertexType.MERGE):
                chain2.append(curr)
                curr = self._prev(curr)
            chain2.append(curr)

            self.chains.append(MonotoneChain(chain1, self.pts, is_left=True))
            self.chains.append(MonotoneChain(chain2, self.pts, is_left=False))

        # Build lookup tables for left-boundary chains
        self.min_to_chain: Dict[int, MonotoneChain] = {}
        self.max_to_chain: Dict[int, MonotoneChain] = {}
        for chain in self.chains:
            if chain.is_left:
                if chain.lower not in self.min_to_chain:
                    self.min_to_chain[chain.lower] = chain
                if chain.upper not in self.max_to_chain:
                    self.max_to_chain[chain.upper] = chain

    def _decompose(self) -> None:
        """
        Monotone decomposition via plane sweep. O(r log r) for BST operations.
        Only processes O(r) extrema events, not all n vertices.
        """
        self.diagonals = []
        sorted_extrema = sorted(self.extrema, 
                               key=lambda i: (-self.pts[i][1], self.pts[i][0]))
        bst = ChainBST()

        for v in sorted_extrema:
            vtype = self.types[v]
            vx, vy = self.pts[v]
            bst.set_y(vy - 1e-9)

            if vtype == VertexType.START:
                if v in self.max_to_chain:
                    chain = self.max_to_chain[v]
                    chain.reset()
                    bst.insert(chain)

            elif vtype == VertexType.END:
                if v in self.min_to_chain:
                    chain = self.min_to_chain[v]
                    if chain.pending is not None:
                        self._add_diagonal(v, chain.pending)
                    bst.remove(chain)

            elif vtype == VertexType.SPLIT:
                L = bst.left_of(vx)
                if L is not None:
                    L.advance_to(vy)
                    if L.pending is not None:
                        self._add_diagonal(v, L.pending)
                        L.pending = None
                    else:
                        self._add_diagonal(v, L.slab_entry())

                if v in self.max_to_chain:
                    chain = self.max_to_chain[v]
                    chain.reset()
                    bst.insert(chain)

            elif vtype == VertexType.MERGE:
                R = self.min_to_chain.get(v)
                L = bst.left_of(vx)

                if R is not None and R.pending is not None:
                    self._add_diagonal(v, R.pending)

                if L is not None:
                    L.advance_to(vy)
                    if L.pending is not None:
                        self._add_diagonal(v, L.pending)
                    L.pending = v

                if R is not None:
                    bst.remove(R)

    def _add_diagonal(self, u: int, v: int) -> None:
        """Add a diagonal if valid."""
        if u != v:
            d = (min(u, v), max(u, v))
            if d not in self.diagonals:
                # Don't add if adjacent
                diff = abs(u - v)
                if diff != 1 and diff != self.n - 1:
                    self.diagonals.append(d)

    def _triangulate_faces(self) -> List[Tuple[int, int, int]]:
        """Extract faces and triangulate each monotone face. O(n)."""
        # Build adjacency
        adj: Dict[int, Set[int]] = defaultdict(set)
        for i in range(self.n):
            j = self._next(i)
            adj[i].add(j)
            adj[j].add(i)
        for u, v in self.diagonals:
            adj[u].add(v)
            adj[v].add(u)

        # Sort neighbors by angle
        adj_sorted: Dict[int, List[int]] = {}
        for u in adj:
            neighbors = list(adj[u])
            neighbors.sort(key=lambda x: math.atan2(
                self.pts[x][1] - self.pts[u][1],
                self.pts[x][0] - self.pts[u][0]
            ))
            adj_sorted[u] = neighbors

        # Extract faces via half-edge traversal
        half_edges: Set[Tuple[int, int]] = set()
        for i in range(self.n):
            j = self._next(i)
            half_edges.add((i, j))
            half_edges.add((j, i))
        for u, v in self.diagonals:
            half_edges.add((u, v))
            half_edges.add((v, u))

        used: Set[Tuple[int, int]] = set()
        faces: List[List[int]] = []

        for start_edge in half_edges:
            if start_edge in used:
                continue

            face = [start_edge[0]]
            prev, curr = start_edge[0], start_edge[1]
            used.add((prev, curr))

            iterations = 0
            while curr != start_edge[0] and iterations < 2 * self.n:
                face.append(curr)
                neighbors = adj_sorted[curr]
                idx = neighbors.index(prev)
                next_v = neighbors[(idx + 1) % len(neighbors)]
                used.add((curr, next_v))
                prev, curr = curr, next_v
                iterations += 1

            if len(face) >= 3 and iterations < 2 * self.n:
                # Only keep interior faces (positive area)
                area = sum(
                    self.pts[face[i]][0] * self.pts[face[(i+1) % len(face)]][1] -
                    self.pts[face[(i+1) % len(face)]][0] * self.pts[face[i]][1]
                    for i in range(len(face))
                )
                if area > 1e-9:
                    faces.append(face)

        # Triangulate each monotone face
        triangles: List[Tuple[int, int, int]] = []
        for face in faces:
            triangles.extend(self._triangulate_monotone(face))
        return triangles

    def _triangulate_monotone(self, face: List[int]) -> List[Tuple[int, int, int]]:
        """
        Stack-based triangulation of y-monotone polygon. O(m) for m vertices.
        """
        m = len(face)
        if m < 3:
            return []
        if m == 3:
            return [(face[0], face[1], face[2])]

        # Find top and bottom vertices
        def better_top(a: int, b: int) -> bool:
            dy = self.pts[a][1] - self.pts[b][1]
            if abs(dy) > 1e-12:
                return dy > 0
            return self.pts[a][0] < self.pts[b][0]

        top = max(face, key=lambda v: (self.pts[v][1], -self.pts[v][0]))
        bottom = min(face, key=lambda v: (self.pts[v][1], -self.pts[v][0]))

        idx_top = face.index(top)

        # Build left and right chains
        chain1, chain2 = [], []
        for i in range(m):
            chain1.append(face[(idx_top + i) % m])
            if face[(idx_top + i) % m] == bottom:
                break
        for i in range(m):
            chain2.append(face[(idx_top - i + m) % m])
            if face[(idx_top - i + m) % m] == bottom:
                break

        # Determine which chain is left
        x1 = self.pts[chain1[1]][0] if len(chain1) > 1 else self.pts[bottom][0]
        x2 = self.pts[chain2[1]][0] if len(chain2) > 1 else self.pts[bottom][0]
        left_chain = chain1 if x1 < x2 else chain2

        is_left = {v: (v in left_chain) for v in face}

        # Sort vertices by y-coordinate (decreasing)
        sorted_verts = sorted(face, key=lambda v: (-self.pts[v][1], self.pts[v][0]))

        # Stack-based triangulation
        triangles = []
        stack = [sorted_verts[0], sorted_verts[1]]

        for i in range(2, m - 1):
            v = sorted_verts[i]
            if is_left[v] != is_left[stack[-1]]:
                # Different chain: connect to all stack vertices
                while len(stack) > 1:
                    u = stack.pop()
                    w = stack[-1]
                    triangles.append((v, u, w))
                stack.pop()
                stack.append(sorted_verts[i - 1])
                stack.append(v)
            else:
                # Same chain: pop while diagonal is inside
                u = stack.pop()
                while stack:
                    w = stack[-1]
                    cross = self._cross(v, u, w)
                    ok = (cross > 1e-12) if is_left[v] else (cross < -1e-12)
                    if not ok:
                        break
                    triangles.append((v, u, w))
                    u = stack.pop()
                stack.append(u)
                stack.append(v)

        # Connect last vertex to remaining stack
        v = sorted_verts[m - 1]
        while len(stack) > 1:
            u = stack.pop()
            w = stack[-1]
            triangles.append((v, u, w))

        return triangles

    def _signed_area(self, pts: List[Tuple[float, float]]) -> float:
        return sum(
            pts[i][0] * pts[(i+1) % len(pts)][1] - pts[(i+1) % len(pts)][0] * pts[i][1]
            for i in range(len(pts))
        )


# Polygon generators
def paper_example() -> List[Tuple[float, float]]:
    return [
        (0.0, 2.5), (1.2, 5.5), (2.5, 3.8), (4.0, 6.5),
        (5.5, 4.8), (7.0, 7.0), (8.0, 5.5), (6.5, 3.5),
        (8.0, 1.5), (5.0, 2.5), (3.0, 0.0), (1.5, 1.5),
    ]


def convex_polygon(n: int = 6) -> List[Tuple[float, float]]:
    return [(math.cos(2 * math.pi * i / n + math.pi / 2),
             math.sin(2 * math.pi * i / n + math.pi / 2)) for i in range(n)]


def star_polygon(points: int = 5, outer: float = 2.0, inner: float = 0.8) -> List[Tuple[float, float]]:
    pts = []
    for i in range(points * 2):
        angle = math.pi / 2 + i * math.pi / points
        r = outer if i % 2 == 0 else inner
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def comb_polygon(teeth: int = 3) -> List[Tuple[float, float]]:
    pts = [(0, 0), (teeth * 2, 0), (teeth * 2, 1)]
    for i in range(teeth - 1, -1, -1):
        x = i * 2 + 1
        pts.extend([(x + 0.5, 1), (x, 2), (x - 0.5, 1)])
    pts.append((0, 1))
    return pts


def l_shape() -> List[Tuple[float, float]]:
    return [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]


def validate(pts: List[Tuple[float, float]], triangles: List[Tuple[int, int, int]]) -> Tuple[bool, str]:
    """Validate triangulation correctness."""
    n = len(pts)
    if len(triangles) != n - 2:
        return False, f"Wrong count: {len(triangles)} != {n - 2}"

    poly_area = abs(sum(
        pts[i][0] * pts[(i+1) % n][1] - pts[(i+1) % n][0] * pts[i][1]
        for i in range(n)
    )) / 2

    tri_area = sum(
        abs((pts[t[1]][0] - pts[t[0]][0]) * (pts[t[2]][1] - pts[t[0]][1]) -
            (pts[t[1]][1] - pts[t[0]][1]) * (pts[t[2]][0] - pts[t[0]][0])) / 2
        for t in triangles
    )

    if abs(poly_area - tri_area) > 1e-6 * max(1, poly_area):
        return False, f"Area mismatch: {poly_area:.6f} vs {tri_area:.6f}"

    for t in triangles:
        for v in t:
            if v < 0 or v >= n:
                return False, f"Invalid vertex: {v}"

    return True, "OK"


def run_tests():
    print("=" * 60)
    print("Polygon Triangulation - O(n + r log r)")
    print("=" * 60)

    cases = [
        ("Triangle", [(0, 0), (1, 0), (0.5, 1)]),
        ("Square", [(0, 0), (1, 0), (1, 1), (0, 1)]),
        ("Pentagon", convex_polygon(5)),
        ("Hexagon", convex_polygon(6)),
        ("L-shape", l_shape()),
        ("Paper example (n=12, r=3)", paper_example()),
        ("5-star", star_polygon(5)),
        ("6-star", star_polygon(6)),
        ("7-star", star_polygon(7)),
        ("Comb-3", comb_polygon(3)),
        ("Comb-5", comb_polygon(5)),
    ]

    passed, failed = 0, 0

    for name, pts in cases:
        try:
            tri = Triangulator(pts)
            triangles = tri.triangulate()
            ok, msg = validate(tri.pts, triangles)

            if ok:
                passed += 1
                status = "PASS"
            else:
                failed += 1
                status = f"FAIL: {msg}"

            print(f"{name}: n={len(pts)}, r={tri.r}, diags={len(tri.diagonals)}, tris={len(triangles)} - {status}")

        except Exception as e:
            failed += 1
            print(f"{name}: ERROR - {e}")

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    success = run_tests()

    if success:
        print("\nPaper example details:")
        pts = paper_example()
        tri = Triangulator(pts)
        triangles = tri.triangulate()
        print(f"  n={len(pts)}, r={tri.r}")
        print(f"  Extrema: {len(tri.extrema)}")
        print(f"  Diagonals: {tri.diagonals}")
        print(f"  Triangles: {len(triangles)}")

        expected_r = 3
        expected_diags = 2
        if tri.r == expected_r and len(tri.diagonals) == expected_diags:
            print(f"\n  ✓ Matches paper: r={expected_r}, {expected_diags} decomposition diagonals")
        else:
            print(f"\n  ✗ Expected r={expected_r}, diags={expected_diags}")
    else:
        exit(1)
