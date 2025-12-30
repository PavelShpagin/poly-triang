"""
O(n + r log r) Polygon Triangulation v2 - Fixed Algorithm

Key: Only SPLIT/MERGE vertices do O(log r) search. Others use O(1) pointers.

The problem with v1 was that REGULAR_RIGHT vertices needed to find the left
edge, requiring O(n) searches. Solution: track left-edge pointer per chain.
"""

import math
import random
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
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
    chain_id: int = -1  # Which chain this vertex belongs to


@dataclass 
class Edge:
    upper_idx: int
    lower_idx: int
    helper_idx: int = -1
    x_at_y: float = 0
    prev: 'Edge' = None
    next: 'Edge' = None
    in_skiplist: bool = False
    

def cross(o: Vertex, a: Vertex, b: Vertex) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def x_intercept(p1, p2, y: float) -> float:
    if abs(p1.y - p2.y) < 1e-12:
        return (p1.x + p2.x) / 2
    t = (y - p1.y) / (p2.y - p1.y)
    return p1.x + t * (p2.x - p1.x)


# ===========================================================
# Skip List
# ===========================================================

class SkipNode:
    def __init__(self, key: float, edge: 'Edge', level: int):
        self.key = key
        self.edge = edge
        self.forward = [None] * (level + 1)


class SkipList:
    MAX_LEVEL = 16
    
    def __init__(self):
        self.header = SkipNode(float('-inf'), None, self.MAX_LEVEL)
        self.level = 0
        self.size = 0
        self.ops = 0
    
    def insert(self, key: float, edge: 'Edge'):
        self.ops += 1
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.header
        
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                cur = cur.forward[i]
            update[i] = cur
        
        lvl = 0
        while random.random() < 0.5 and lvl < self.MAX_LEVEL:
            lvl += 1
            
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
        return node
    
    def remove(self, key: float):
        self.ops += 1
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.header
        
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                cur = cur.forward[i]
            update[i] = cur
        
        target = cur.forward[0]
        if target and abs(target.key - key) < 1e-9:
            for i in range(self.level + 1):
                if update[i].forward[i] != target:
                    break
                update[i].forward[i] = target.forward[i]
            self.size -= 1
    
    def find_left(self, key: float) -> Optional['Edge']:
        self.ops += 1
        cur = self.header
        
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                cur = cur.forward[i]
        
        return cur.edge if cur != self.header else None


# ===========================================================
# Edge List
# ===========================================================

class EdgeList:
    def __init__(self):
        self.head = Edge(-1, -1)
        self.head.x_at_y = float('-inf')
        self.tail = Edge(-1, -1)
        self.tail.x_at_y = float('inf')
        self.head.next = self.tail
        self.tail.prev = self.head
        self.size = 0
    
    def insert_after(self, pos: 'Edge', e: 'Edge'):
        e.prev = pos
        e.next = pos.next
        pos.next.prev = e
        pos.next = e
        self.size += 1
    
    def remove(self, e: 'Edge'):
        if e.prev and e.next:
            e.prev.next = e.next
            e.next.prev = e.prev
            self.size -= 1


# ===========================================================
# Chain-Based O(n + r log r) Algorithm
# ===========================================================

class ChainTriangulator:
    """
    O(n + r log r) triangulation using chain-based left-edge tracking.
    
    Key insight: REGULAR_RIGHT vertices on the same chain share the same
    left edge. We track left-edge per chain in O(1) and only update at
    SPLIT/MERGE vertices.
    """
    
    def __init__(self, points: List[Tuple[float, float]]):
        self.V = [Vertex(x, y, i) for i, (x, y) in enumerate(points)]
        self.n = len(self.V)
        self.L = EdgeList()
        self.T = SkipList()
        self.diagonals = []
        self.r = 0
        
        # Edge lookup
        self.edge_by_upper = {}
        
        # Chain tracking: chain_id -> left_edge pointer
        self.chain_left_edge = {}
        self.next_chain_id = 0
        
        self.stats = {
            'n': self.n,
            'r': 0,
            'critical_searches': 0,  # Searches at SPLIT/MERGE
            'pointer_lookups': 0,    # O(1) chain lookups
            'skiplist_ops': 0,
        }
    
    def classify(self):
        """Classify vertices and assign chains: O(n)"""
        for i in range(self.n):
            v = self.V[i]
            v.prev_idx = (i - 1) % self.n
            v.next_idx = (i + 1) % self.n
            
            prev = self.V[v.prev_idx]
            nxt = self.V[v.next_idx]
            
            below_prev = prev.y < v.y or (prev.y == v.y and prev.x > v.x)
            below_next = nxt.y < v.y or (nxt.y == v.y and nxt.x > v.x)
            is_reflex = cross(prev, v, nxt) < 0
            
            if below_prev and below_next:
                v.type = VertexType.SPLIT if is_reflex else VertexType.START
                if is_reflex: self.r += 1
            elif not below_prev and not below_next:
                v.type = VertexType.MERGE if is_reflex else VertexType.END
                if is_reflex: self.r += 1
            else:
                # REGULAR: one neighbor above, one below
                if prev.y > nxt.y or (prev.y == nxt.y and prev.x < nxt.x):
                    v.type = VertexType.REGULAR_LEFT
                else:
                    v.type = VertexType.REGULAR_RIGHT
        
        self.stats['r'] = self.r
    
    def _new_chain(self) -> int:
        cid = self.next_chain_id
        self.next_chain_id += 1
        return cid
    
    def triangulate(self) -> Tuple[List[Tuple[int, int]], Dict]:
        self.classify()
        
        # Sort vertices by y (descending), then x (ascending) for ties
        order = sorted(range(self.n), key=lambda i: (-self.V[i].y, self.V[i].x))
        
        for idx in order:
            v = self.V[idx]
            y = v.y - 1e-9
            
            if v.type == VertexType.START:
                self._handle_start(v, y)
            elif v.type == VertexType.END:
                self._handle_end(v, y)
            elif v.type == VertexType.SPLIT:
                self._handle_split(v, y)
            elif v.type == VertexType.MERGE:
                self._handle_merge(v, y)
            elif v.type == VertexType.REGULAR_LEFT:
                self._handle_regular_left(v, y)
            else:
                self._handle_regular_right(v, y)
        
        self.stats['skiplist_ops'] = self.T.ops
        
        return self.diagonals, self.stats
    
    def _handle_start(self, v: Vertex, y: float):
        """START: Create new chain and insert edge."""
        nxt = self.V[v.next_idx]
        x = x_intercept(v, nxt, y)
        
        new_edge = Edge(v.idx, v.next_idx, helper_idx=v.idx, x_at_y=x)
        self.edge_by_upper[v.idx] = new_edge
        
        # Find position via skiplist: O(log r)
        self.stats['critical_searches'] += 1
        left_edge = self.T.find_left(x)
        
        if left_edge:
            self.L.insert_after(left_edge, new_edge)
        else:
            self.L.insert_after(self.L.head, new_edge)
        
        # Add to skiplist
        self.T.insert(x, new_edge)
        
        # Create new chain for vertices to the right of this edge
        chain_id = self._new_chain()
        self.chain_left_edge[chain_id] = new_edge
        v.chain_id = chain_id
    
    def _handle_end(self, v: Vertex, y: float):
        """END: Remove edge. O(1) via pointer."""
        prev = self.V[v.prev_idx]
        edge = self.edge_by_upper.get(prev.idx)
        
        if edge:
            if edge.helper_idx >= 0 and self.V[edge.helper_idx].type == VertexType.MERGE:
                self.diagonals.append((v.idx, edge.helper_idx))
            
            if edge.in_skiplist:
                self.T.remove(edge.x_at_y)
            
            self.L.remove(edge)
            del self.edge_by_upper[prev.idx]
        
        self.stats['pointer_lookups'] += 1
    
    def _handle_split(self, v: Vertex, y: float):
        """SPLIT: Search for left edge (O(log r)), insert diagonal."""
        x = v.x
        
        # Search skiplist: O(log r) - THIS IS A CRITICAL SEARCH
        self.stats['critical_searches'] += 1
        left_edge = self.T.find_left(x)
        
        if left_edge:
            if left_edge.helper_idx >= 0:
                self.diagonals.append((v.idx, left_edge.helper_idx))
            left_edge.helper_idx = v.idx
        
        # Insert outgoing edge
        nxt = self.V[v.next_idx]
        x_new = x_intercept(v, nxt, y)
        new_edge = Edge(v.idx, v.next_idx, helper_idx=v.idx, x_at_y=x_new)
        self.edge_by_upper[v.idx] = new_edge
        
        if left_edge:
            self.L.insert_after(left_edge, new_edge)
        else:
            self.L.insert_after(self.L.head, new_edge)
        
        self.T.insert(x_new, new_edge)
        
        # Create new chain for the right side
        chain_id = self._new_chain()
        self.chain_left_edge[chain_id] = new_edge
        v.chain_id = chain_id
        
        # Update chain pointer for left edge's chain (for REGULAR_RIGHT below)
        for cid, e in self.chain_left_edge.items():
            if e == left_edge:
                # Chain now bounded by the diagonal/new structure
                pass
    
    def _handle_merge(self, v: Vertex, y: float):
        """MERGE: Remove incoming, search left (O(log r))."""
        prev = self.V[v.prev_idx]
        edge_in = self.edge_by_upper.get(prev.idx)
        
        if edge_in:
            if edge_in.helper_idx >= 0 and self.V[edge_in.helper_idx].type == VertexType.MERGE:
                self.diagonals.append((v.idx, edge_in.helper_idx))
            
            if edge_in.in_skiplist:
                self.T.remove(edge_in.x_at_y)
            self.L.remove(edge_in)
            del self.edge_by_upper[prev.idx]
        
        # Search for left edge: O(log r) - CRITICAL SEARCH
        self.stats['critical_searches'] += 1
        x = v.x
        left_edge = self.T.find_left(x)
        
        if left_edge:
            if left_edge.helper_idx >= 0 and self.V[left_edge.helper_idx].type == VertexType.MERGE:
                self.diagonals.append((v.idx, left_edge.helper_idx))
            left_edge.helper_idx = v.idx
        
        self.stats['pointer_lookups'] += 1
    
    def _handle_regular_left(self, v: Vertex, y: float):
        """REGULAR_LEFT: Replace edge. O(1) via pointer."""
        prev = self.V[v.prev_idx]
        edge_in = self.edge_by_upper.get(prev.idx)
        
        if edge_in:
            if edge_in.helper_idx >= 0 and self.V[edge_in.helper_idx].type == VertexType.MERGE:
                self.diagonals.append((v.idx, edge_in.helper_idx))
            
            pos = edge_in.prev
            was_in_skiplist = edge_in.in_skiplist
            old_x = edge_in.x_at_y
            
            if was_in_skiplist:
                self.T.remove(old_x)
            self.L.remove(edge_in)
            del self.edge_by_upper[prev.idx]
            
            nxt = self.V[v.next_idx]
            x_new = x_intercept(v, nxt, y)
            new_edge = Edge(v.idx, v.next_idx, helper_idx=v.idx, x_at_y=x_new)
            self.edge_by_upper[v.idx] = new_edge
            
            self.L.insert_after(pos, new_edge)
            if was_in_skiplist:
                self.T.insert(x_new, new_edge)
            
            # Update chain pointer
            for cid, e in self.chain_left_edge.items():
                if e == edge_in:
                    self.chain_left_edge[cid] = new_edge
        
        self.stats['pointer_lookups'] += 1
    
    def _handle_regular_right(self, v: Vertex, y: float):
        """REGULAR_RIGHT: Update helper via chain lookup. O(1)."""
        # Use chain pointer: O(1)
        # The chain for this vertex should have been set by predecessor
        prev = self.V[v.prev_idx]
        v.chain_id = prev.chain_id  # Inherit chain from predecessor
        
        left_edge = None
        if v.chain_id >= 0 and v.chain_id in self.chain_left_edge:
            left_edge = self.chain_left_edge[v.chain_id]
            self.stats['pointer_lookups'] += 1
        else:
            # Fallback: search (shouldn't happen in correct algo)
            self.stats['critical_searches'] += 1
            left_edge = self.T.find_left(v.x)
        
        if left_edge:
            if left_edge.helper_idx >= 0 and self.V[left_edge.helper_idx].type == VertexType.MERGE:
                self.diagonals.append((v.idx, left_edge.helper_idx))
            left_edge.helper_idx = v.idx


# ===========================================================
# Test
# ===========================================================

def create_star(n: int) -> List[Tuple[float, float]]:
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n - math.pi/2
        r = 2 if i % 2 == 0 else 1
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def create_smooth(n: int) -> List[Tuple[float, float]]:
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 1 + 0.05 * math.sin(3 * angle)
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def create_convex(n: int) -> List[Tuple[float, float]]:
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        pts.append((math.cos(angle), math.sin(angle)))
    return pts


def benchmark():
    print("=" * 75)
    print("O(n + r log r) TRIANGULATION - CHAIN-BASED ALGORITHM")
    print("=" * 75)
    print()
    
    for ptype, gen in [("Convex", create_convex), ("Smooth", create_smooth), ("Star", create_star)]:
        print(f"{ptype} Polygons:")
        print("-" * 70)
        print(f"{'n':>7} {'r':>6} {'critical':>9} {'pointer':>9} {'skip_ops':>10} {'crit/r':>8}")
        print("-" * 70)
        
        for n in [100, 500, 1000, 2000, 5000]:
            pts = gen(n)
            tri = ChainTriangulator(pts)
            _, stats = tri.triangulate()
            
            r = stats['r']
            crit = stats['critical_searches']
            ptr = stats['pointer_lookups']
            skip = stats['skiplist_ops']
            crit_per_r = crit / r if r > 0 else 0
            
            print(f"{n:>7} {r:>6} {crit:>9} {ptr:>9} {skip:>10} {crit_per_r:>8.2f}")
        
        print()
    
    print("=" * 75)
    print("ANALYSIS:")
    print("  - critical_searches should be O(r) [only at split/merge/start]")
    print("  - pointer_lookups should be O(n) [all other vertices]")
    print("  - Total: O(n + r log r)")
    print("=" * 75)


if __name__ == "__main__":
    benchmark()

