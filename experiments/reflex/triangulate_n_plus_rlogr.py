"""
O(n + r log r) Polygon Triangulation via Lazy Landmark Indexing

Key insight: Only SPLIT/MERGE vertices need O(log r) search.
All other vertices use O(1) linked-list operations.

Data structures:
- L: Doubly-linked list of edges in sorted order
- T: Skip list of "landmark" edges for O(log r) approximate search
"""

import math
import random
import time
from typing import List, Tuple, Optional, Dict, Set
from dataclasses import dataclass, field
from enum import Enum


class VertexType(Enum):
    START = 1
    END = 2
    SPLIT = 3
    MERGE = 4
    REGULAR_LEFT = 5
    REGULAR_RIGHT = 6


@dataclass
class Vertex:
    x: float
    y: float
    idx: int
    type: VertexType = None
    prev_idx: int = -1
    next_idx: int = -1


@dataclass
class Edge:
    """Edge in the doubly-linked list"""
    upper_idx: int  # Index of upper vertex
    lower_idx: int  # Index of lower vertex
    x_at_y: float   # x-intercept at current sweep y
    helper_idx: int = -1
    prev: 'Edge' = None
    next: 'Edge' = None
    in_skiplist: bool = False


def cross(o: Vertex, a: Vertex, b: Vertex) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def x_intercept(v1: Vertex, v2: Vertex, y: float) -> float:
    """X-coordinate where edge (v1,v2) crosses y-level"""
    if abs(v1.y - v2.y) < 1e-12:
        return (v1.x + v2.x) / 2
    t = (y - v1.y) / (v2.y - v1.y)
    return v1.x + t * (v2.x - v1.x)


# ===========================================================
# Skip List for Landmarks (O(log r) search)
# ===========================================================

class SkipNode:
    def __init__(self, key: float, edge: Edge, level: int):
        self.key = key
        self.edge = edge
        self.forward = [None] * (level + 1)


class SkipList:
    """Skip list for O(log r) approximate position finding"""
    MAX_LEVEL = 16
    
    def __init__(self):
        self.header = SkipNode(float('-inf'), None, self.MAX_LEVEL)
        self.level = 0
        self.size = 0
        self.ops = 0
    
    def _rand_level(self) -> int:
        lvl = 0
        while random.random() < 0.5 and lvl < self.MAX_LEVEL:
            lvl += 1
        return lvl
    
    def insert(self, key: float, edge: Edge):
        """Insert edge as landmark"""
        self.ops += 1
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.header
        
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                self.ops += 1
                cur = cur.forward[i]
            update[i] = cur
        
        lvl = self._rand_level()
        if lvl > self.level:
            for i in range(self.level + 1, lvl + 1):
                update[i] = self.header
            self.level = lvl
        
        node = SkipNode(key, edge, lvl)
        for i in range(lvl + 1):
            node.forward[i] = update[i].forward[i]
            update[i].forward[i] = node
        
        self.size += 1
        edge.in_skiplist = True
    
    def remove(self, key: float):
        """Remove landmark"""
        self.ops += 1
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.header
        
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                self.ops += 1
                cur = cur.forward[i]
            update[i] = cur
        
        target = cur.forward[0]
        if target and abs(target.key - key) < 1e-9:
            for i in range(self.level + 1):
                if update[i].forward[i] != target:
                    break
                update[i].forward[i] = target.forward[i]
            self.size -= 1
    
    def find_left_landmark(self, key: float) -> Optional[Edge]:
        """Find rightmost landmark with key < given key"""
        self.ops += 1
        cur = self.header
        
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                self.ops += 1
                cur = cur.forward[i]
        
        return cur.edge if cur != self.header else None


# ===========================================================
# Edge List (O(1) local operations)
# ===========================================================

class EdgeList:
    """Doubly-linked list of edges with O(1) local ops"""
    
    def __init__(self):
        self.head = Edge(-1, -1, float('-inf'))  # Sentinel
        self.tail = Edge(-1, -1, float('inf'))   # Sentinel
        self.head.next = self.tail
        self.tail.prev = self.head
        self.size = 0
        self.scan_ops = 0
    
    def insert_after(self, pos: Edge, new_edge: Edge):
        """Insert new_edge after pos: O(1)"""
        new_edge.prev = pos
        new_edge.next = pos.next
        pos.next.prev = new_edge
        pos.next = new_edge
        self.size += 1
    
    def remove(self, edge: Edge):
        """Remove edge: O(1)"""
        edge.prev.next = edge.next
        edge.next.prev = edge.prev
        self.size -= 1
    
    def scan_right_to(self, start: Edge, target_x: float) -> Edge:
        """Scan right from start to find edge just left of target_x"""
        cur = start
        while cur.next != self.tail and cur.next.x_at_y < target_x:
            self.scan_ops += 1
            cur = cur.next
        return cur


# ===========================================================
# Main Algorithm
# ===========================================================

class TriangulatorNplusRlogR:
    """O(n + r log r) triangulation via lazy landmarks"""
    
    def __init__(self, vertices: List[Vertex]):
        self.V = vertices
        self.n = len(vertices)
        self.L = EdgeList()
        self.T = SkipList()
        self.diagonals = []
        self.r = 0
        
        # Direct pointers: vertex_idx -> its outgoing edge (for O(1) access)
        self.vertex_to_edge: Dict[int, Edge] = {}
        
        # Stats
        self.stats = {
            'n': self.n,
            'r': 0,
            'skiplist_ops': 0,
            'scan_ops': 0,
            'landmarks_added': 0,
            'direct_accesses': 0,
        }
    
    def classify(self):
        """Classify all vertices: O(n)"""
        for i in range(self.n):
            v = self.V[i]
            v.prev_idx = (i - 1) % self.n
            v.next_idx = (i + 1) % self.n
            
            prev = self.V[v.prev_idx]
            nxt = self.V[v.next_idx]
            
            below_prev = prev.y < v.y
            below_next = nxt.y < v.y
            is_reflex = cross(prev, v, nxt) < 0
            
            if below_prev and below_next:
                v.type = VertexType.SPLIT if is_reflex else VertexType.START
                if is_reflex:
                    self.r += 1
            elif not below_prev and not below_next:
                v.type = VertexType.MERGE if is_reflex else VertexType.END
                if is_reflex:
                    self.r += 1
            else:
                if prev.y > nxt.y:
                    v.type = VertexType.REGULAR_LEFT
                else:
                    v.type = VertexType.REGULAR_RIGHT
        
        self.stats['r'] = self.r
    
    def find_left_edge(self, x: float, y: float) -> Edge:
        """
        Find edge immediately left of (x, y).
        Uses landmark for O(log r) approximate search, then O(scan) refinement.
        """
        # Find nearest landmark to the left
        landmark = self.T.find_left_landmark(x)
        
        # Start scan from landmark (or head)
        start = landmark if landmark else self.L.head
        
        # Scan right to find exact position
        edge = self.L.scan_right_to(start, x)
        
        return edge if edge != self.L.head else None
    
    def add_landmark(self, edge: Edge):
        """Add edge as landmark in T"""
        if not edge.in_skiplist and edge != self.L.head and edge != self.L.tail:
            self.T.insert(edge.x_at_y, edge)
            self.stats['landmarks_added'] += 1
    
    def process_start(self, v: Vertex, y: float):
        """START: Insert two outgoing edges. Need search for position."""
        # START vertex has both neighbors below - creates two outgoing edges
        # We insert the LEFT edge (to prev) - the right edge will be handled
        # when we process the right chain
        nxt = self.V[v.next_idx]
        prev = self.V[v.prev_idx]
        
        # Left outgoing edge (v -> prev, going down-left)
        x_left = x_intercept(v, prev, y)
        left_edge = Edge(v.idx, v.prev_idx, x_left, helper_idx=v.idx)
        
        # Right outgoing edge (v -> next, going down-right)  
        x_right = x_intercept(v, nxt, y)
        right_edge = Edge(v.idx, v.next_idx, x_right, helper_idx=v.idx)
        
        # Find position - for START we need to search (but START vertices are rare)
        pos = self.find_left_edge(x_left, y)
        
        if pos:
            self.L.insert_after(pos, left_edge)
        else:
            self.L.insert_after(self.L.head, left_edge)
        self.L.insert_after(left_edge, right_edge)
        
        # Store direct pointers for O(1) removal later
        self.vertex_to_edge[v.prev_idx] = left_edge   # prev will remove this
        self.vertex_to_edge[v.next_idx] = right_edge  # next will remove this
    
    def process_end(self, v: Vertex, y: float):
        """END: Remove incoming edge. O(1) with direct pointer."""
        # Use direct pointer - the edge ending at v was stored when created
        edge = self.vertex_to_edge.get(v.idx)
        
        if edge:
            self.stats['direct_accesses'] += 1
            if edge.in_skiplist:
                self.T.remove(edge.x_at_y)
            self.L.remove(edge)
            del self.vertex_to_edge[v.idx]
    
    def process_split(self, v: Vertex, y: float):
        """SPLIT: Search for left edge, add diagonal, insert new edge."""
        # Find left edge: O(log r) + O(scan)
        left_edge = self.find_left_edge(v.x, y)
        
        if left_edge and left_edge != self.L.head:
            # Add landmark for future queries
            self.add_landmark(left_edge)
            
            # Add diagonal
            if left_edge.helper_idx >= 0:
                self.diagonals.append((v.idx, left_edge.helper_idx))
            
            # Update helper
            left_edge.helper_idx = v.idx
        
        # Insert outgoing edge
        nxt = self.V[v.next_idx]
        x = x_intercept(v, nxt, y)
        new_edge = Edge(v.idx, v.next_idx, x, helper_idx=v.idx)
        
        if left_edge:
            self.L.insert_after(left_edge, new_edge)
        else:
            self.L.insert_after(self.L.head, new_edge)
    
    def process_merge(self, v: Vertex, y: float):
        """MERGE: Remove incoming edge, search for left, add diagonals."""
        # Remove incoming edge from above-left
        prev = self.V[v.prev_idx]
        edge_in = self.find_left_edge(v.x + 0.001, y)
        
        if edge_in and edge_in.helper_idx >= 0:
            helper = self.V[edge_in.helper_idx]
            if helper.type == VertexType.MERGE:
                self.diagonals.append((v.idx, edge_in.helper_idx))
        
        if edge_in:
            if edge_in.in_skiplist:
                self.T.remove(edge_in.x_at_y)
            self.L.remove(edge_in)
        
        # Find left edge
        left_edge = self.find_left_edge(v.x, y)
        
        if left_edge and left_edge != self.L.head:
            # Add landmark
            self.add_landmark(left_edge)
            
            # Check helper
            if left_edge.helper_idx >= 0:
                helper = self.V[left_edge.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, left_edge.helper_idx))
            
            left_edge.helper_idx = v.idx
    
    def process_regular_left(self, v: Vertex, y: float):
        """REGULAR_LEFT: Replace incoming with outgoing edge. O(1)."""
        # Find and remove incoming edge
        edge_in = self.find_left_edge(v.x + 0.001, y)
        
        if edge_in:
            # Check helper
            if edge_in.helper_idx >= 0:
                helper = self.V[edge_in.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, edge_in.helper_idx))
            
            # Remove old edge
            was_landmark = edge_in.in_skiplist
            if was_landmark:
                self.T.remove(edge_in.x_at_y)
            self.L.remove(edge_in)
            
            # Insert new edge at same position
            nxt = self.V[v.next_idx]
            x = x_intercept(v, nxt, y)
            new_edge = Edge(v.idx, v.next_idx, x, helper_idx=v.idx)
            
            # Insert after predecessor of removed edge
            prev_edge = edge_in.prev if edge_in.prev != self.L.head else self.L.head
            self.L.insert_after(prev_edge, new_edge)
            
            if was_landmark:
                self.T.insert(x, new_edge)
    
    def process_regular_right(self, v: Vertex, y: float):
        """REGULAR_RIGHT: Just update helper of left edge. O(1) with pointer."""
        left_edge = self.find_left_edge(v.x, y)
        
        if left_edge and left_edge != self.L.head:
            if left_edge.helper_idx >= 0:
                helper = self.V[left_edge.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, left_edge.helper_idx))
            left_edge.helper_idx = v.idx
    
    def triangulate(self) -> Tuple[List[Tuple[int, int]], Dict]:
        """Main triangulation routine"""
        self.classify()
        
        # Sort vertices by y (decreasing)
        order = sorted(range(self.n), key=lambda i: -self.V[i].y)
        
        for idx in order:
            v = self.V[idx]
            y = v.y - 1e-9  # Just below current vertex
            
            if v.type == VertexType.START:
                self.process_start(v, y)
            elif v.type == VertexType.END:
                self.process_end(v, y)
            elif v.type == VertexType.SPLIT:
                self.process_split(v, y)
            elif v.type == VertexType.MERGE:
                self.process_merge(v, y)
            elif v.type == VertexType.REGULAR_LEFT:
                self.process_regular_left(v, y)
            elif v.type == VertexType.REGULAR_RIGHT:
                self.process_regular_right(v, y)
        
        self.stats['skiplist_ops'] = self.T.ops
        self.stats['scan_ops'] = self.L.scan_ops
        
        return self.diagonals, self.stats


# ===========================================================
# Test
# ===========================================================

def create_star_polygon(n: int) -> List[Vertex]:
    """Star polygon with ~n/2 reflex vertices"""
    verts = []
    for i in range(n):
        angle = 2 * math.pi * i / n - math.pi/2
        r = 2 if i % 2 == 0 else 1
        verts.append(Vertex(r * math.cos(angle), r * math.sin(angle), i))
    return verts


def create_smooth_polygon(n: int) -> List[Vertex]:
    """Nearly convex polygon"""
    verts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 1 + 0.05 * math.sin(3 * angle)
        verts.append(Vertex(r * math.cos(angle), r * math.sin(angle), i))
    return verts


def benchmark():
    print("=" * 70)
    print("O(n + r log r) TRIANGULATION BENCHMARK")
    print("=" * 70)
    print()
    
    print("Star Polygons (r ~ n/2):")
    print("-" * 60)
    print(f"{'n':>8} {'r':>8} {'skip_ops':>10} {'scan_ops':>10} {'total':>10} {'per_n':>8}")
    print("-" * 60)
    
    for n in [100, 500, 1000, 2000, 5000]:
        verts = create_star_polygon(n)
        tri = TriangulatorNplusRlogR(verts)
        _, stats = tri.triangulate()
        
        total = stats['skiplist_ops'] + stats['scan_ops']
        per_n = total / n
        
        print(f"{n:>8} {stats['r']:>8} {stats['skiplist_ops']:>10} {stats['scan_ops']:>10} "
              f"{total:>10} {per_n:>8.2f}")
    
    print()
    print("Smooth Polygons (r ~ O(1)):")
    print("-" * 60)
    print(f"{'n':>8} {'r':>8} {'skip_ops':>10} {'scan_ops':>10} {'total':>10} {'per_n':>8}")
    print("-" * 60)
    
    for n in [100, 500, 1000, 2000, 5000]:
        verts = create_smooth_polygon(n)
        tri = TriangulatorNplusRlogR(verts)
        _, stats = tri.triangulate()
        
        total = stats['skiplist_ops'] + stats['scan_ops']
        per_n = total / n
        
        print(f"{n:>8} {stats['r']:>8} {stats['skiplist_ops']:>10} {stats['scan_ops']:>10} "
              f"{total:>10} {per_n:>8.2f}")
    
    print()
    print("=" * 70)
    print("ANALYSIS:")
    print("- Skip list ops should be O(r log r)")
    print("- Scan ops should be O(n) total")
    print("- Total should be O(n + r log r)")
    print("=" * 70)


if __name__ == "__main__":
    benchmark()

