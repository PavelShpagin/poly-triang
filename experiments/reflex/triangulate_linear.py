
"""
TRUE LINEAR TIME TRIANGULATION
Combines:
1. Bucket Sort (Radix Sort) for O(n) vertex sorting.
2. Splay Tree Sweep for O(n) triangulation.
"""

import math
import random
import sys
import time
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum

# Increase recursion depth for deep trees
sys.setrecursionlimit(20000)

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
    upper: Vertex
    lower: Vertex
    helper_idx: int
    
    def get_x(self, y: float) -> float:
        if abs(self.upper.y - self.lower.y) < 1e-12:
            return self.upper.x
        t = (y - self.upper.y) / (self.lower.y - self.upper.y)
        return self.upper.x + t * (self.lower.x - self.upper.x)

class SplayNode:
    def __init__(self, edge: Edge):
        self.edge = edge
        self.left = None
        self.right = None
        self.parent = None

class SplayTree:
    def __init__(self):
        self.root = None
        self.splay_ops = 0
        self.comparisons = 0  # Count node visits/comparisons
    
    def _rotate_left(self, x):
        self.splay_ops += 1 # Count rotation as operation unit
        y = x.right
        x.right = y.left
        if y.left:
            y.left.parent = x
        y.parent = x.parent
        if not x.parent:
            self.root = y
        elif x == x.parent.left:
            x.parent.left = y
        else:
            x.parent.right = y
        y.left = x
        x.parent = y

    def _rotate_right(self, x):
        self.splay_ops += 1 # Count rotation as operation unit
        y = x.left
        x.left = y.right
        if y.right:
            y.right.parent = x
        y.parent = x.parent
        if not x.parent:
            self.root = y
        elif x == x.parent.right:
            x.parent.right = y
        else:
            x.parent.left = y
        y.right = x
        x.parent = y

    def splay(self, x):
        if not x: return
        # Do not count splay call itself, count rotations
        while x.parent:
            if not x.parent.parent:
                if x == x.parent.left:
                    self._rotate_right(x.parent)
                else:
                    self._rotate_left(x.parent)
            elif x == x.parent.left and x.parent == x.parent.parent.left:
                self._rotate_right(x.parent.parent)
                self._rotate_right(x.parent)
            elif x == x.parent.right and x.parent == x.parent.parent.right:
                self._rotate_left(x.parent.parent)
                self._rotate_left(x.parent)
            elif x == x.parent.right and x.parent == x.parent.parent.left:
                self._rotate_left(x.parent)
                self._rotate_right(x.parent)
            else:
                self._rotate_right(x.parent)
                self._rotate_left(x.parent)

    def insert(self, edge: Edge, current_y: float):
        if not self.root:
            self.root = SplayNode(edge)
            return

        node = self.root
        x_val = edge.get_x(current_y)
        
        while True:
            self.comparisons += 1
            node_x = node.edge.get_x(current_y)
            if x_val < node_x:
                if node.left:
                    node = node.left
                else:
                    node.left = SplayNode(edge)
                    node.left.parent = node
                    self.splay(node.left)
                    return
            else:
                if node.right:
                    node = node.right
                else:
                    node.right = SplayNode(edge)
                    node.right.parent = node
                    self.splay(node.right)
                    return

    def find_left_edge(self, x: float, current_y: float) -> Optional[Edge]:
        """Find edge immediately to the left of x at height current_y"""
        if not self.root:
            return None
            
        node = self.root
        best = None
        best_node = None
        
        while node:
            self.comparisons += 1
            node_x = node.edge.get_x(current_y)
            if node_x < x:
                best = node.edge
                best_node = node
                node = node.right
            else:
                node = node.left
        
        if best_node:
            self.splay(best_node)
            
        return best

    def delete(self, edge: Edge, current_y: float):
        if not self.root:
            return

        node = self.root
        target_x = edge.get_x(current_y)
        target_node = None
        
        while node:
            self.comparisons += 1
            if node.edge == edge: 
                target_node = node
                break
            node_x = node.edge.get_x(current_y)
            if target_x < node_x:
                node = node.left
            else:
                node = node.right
        
        if not target_node:
            return 

        self.splay(target_node) 
        
        if not target_node.left:
            self.root = target_node.right
            if self.root:
                self.root.parent = None
        else:
            right_tree = target_node.right
            self.root = target_node.left
            self.root.parent = None
            
            m = self.root
            while m.right:
                self.comparisons += 1
                m = m.right
            self.splay(m) 
            m.right = right_tree
            if right_tree:
                right_tree.parent = m

def cross(o, a, b):
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)

class LinearTriangulator:
    def __init__(self, points: List[Tuple[float, float]]):
        self.raw_points = points
        self.n = len(points)
        self.V = []
        self.diagonals = []
        self.tree = SplayTree()
        self.stats = {'n': self.n, 'splay_ops': 0, 'comparisons': 0, 'sort_time': 0, 'sweep_time': 0}
        self.edge_map = {} 
 

    def _bucket_sort_vertices(self):
        """Sort vertices by Y (descending) then X using buckets."""
        if self.n == 0: return []
        
        min_y = min(p[1] for p in self.raw_points)
        max_y = max(p[1] for p in self.raw_points)
        rng = max_y - min_y
        
        if rng < 1e-9: # All same Y?
            return list(range(self.n))
            
        # Create buckets
        # For N points, N buckets is optimal for uniform distribution.
        bucket_count = self.n
        buckets = [[] for _ in range(bucket_count + 1)]
        
        for i in range(self.n):
            y = self.raw_points[i][1]
            # Invert y because we want descending sort
            # But simpler to sort ascending and reverse?
            # Let's map max_y to 0 and min_y to bucket_count
            
            # We want descending Y.
            # norm = (max_y - y) / rng  => 0 to 1
            # idx = int(norm * bucket_count)
            
            val = (max_y - y) / rng
            b_idx = int(val * bucket_count)
            if b_idx > bucket_count: b_idx = bucket_count
            buckets[b_idx].append(i)
            
        # Collect and sort within buckets (Insertion sort or just standard sort if small)
        # For general inputs, Python's Timsort on small lists is very fast (O(k)).
        sorted_indices = []
        for b in buckets:
            if not b: continue
            # Sort by -y, then x
            b.sort(key=lambda i: (-self.raw_points[i][1], self.raw_points[i][0]))
            sorted_indices.extend(b)
            
        return sorted_indices

    def classify(self):
        self.V = [Vertex(x, y, i) for i, (x, y) in enumerate(self.raw_points)]
        for i in range(self.n):
            v = self.V[i]
            v.prev_idx = (i - 1) % self.n
            v.next_idx = (i + 1) % self.n
            
            prev, nxt = self.V[v.prev_idx], self.V[v.next_idx]
            bp = prev.y < v.y or (prev.y == v.y and prev.x > v.x)
            bn = nxt.y < v.y or (nxt.y == v.y and nxt.x > v.x)
            reflex = cross(prev, v, nxt) < 0
            
            if bp and bn:
                v.type = VertexType.SPLIT if reflex else VertexType.START
            elif not bp and not bn:
                v.type = VertexType.MERGE if reflex else VertexType.END
            else:
                v.type = VertexType.REGULAR_LEFT if not bp else VertexType.REGULAR_RIGHT

    def triangulate(self):
        t0 = time.time()
        
        self.classify()
        
        # Linear Sort
        sorted_indices = self._bucket_sort_vertices()
        
        t1 = time.time()
        self.stats['sort_time'] = t1 - t0
        
        # Sweep
        for idx in sorted_indices:
            v = self.V[idx]
            y = v.y
            
            if v.type == VertexType.START:
                nxt = self.V[v.next_idx]
                edge = Edge(v, nxt, v.idx)
                self.edge_map[v.idx] = edge
                self.tree.insert(edge, y)
                
            elif v.type == VertexType.END:
                if v.prev_idx in self.edge_map:
                    edge = self.edge_map[v.prev_idx]
                    if self.V[edge.helper_idx].type == VertexType.MERGE:
                        self.diagonals.append((v.idx, edge.helper_idx))
                    self.tree.delete(edge, y)
                    del self.edge_map[v.prev_idx]
                    
            elif v.type == VertexType.SPLIT:
                left_edge = self.tree.find_left_edge(v.x, y)
                if left_edge:
                    self.diagonals.append((v.idx, left_edge.helper_idx))
                    left_edge.helper_idx = v.idx
                
                nxt = self.V[v.next_idx]
                edge = Edge(v, nxt, v.idx)
                self.edge_map[v.idx] = edge
                self.tree.insert(edge, y)
                
            elif v.type == VertexType.MERGE:
                if v.prev_idx in self.edge_map:
                    edge = self.edge_map[v.prev_idx]
                    if self.V[edge.helper_idx].type == VertexType.MERGE:
                        self.diagonals.append((v.idx, edge.helper_idx))
                    self.tree.delete(edge, y)
                    del self.edge_map[v.prev_idx]
                
                left_edge = self.tree.find_left_edge(v.x, y)
                if left_edge:
                    if self.V[left_edge.helper_idx].type == VertexType.MERGE:
                        self.diagonals.append((v.idx, left_edge.helper_idx))
                    left_edge.helper_idx = v.idx
                    
            elif v.type == VertexType.REGULAR_LEFT:
                if v.prev_idx in self.edge_map:
                    edge = self.edge_map[v.prev_idx]
                    if self.V[edge.helper_idx].type == VertexType.MERGE:
                        self.diagonals.append((v.idx, edge.helper_idx))
                    self.tree.delete(edge, y)
                    del self.edge_map[v.prev_idx]
                    
                nxt = self.V[v.next_idx]
                new_edge = Edge(v, nxt, v.idx)
                self.edge_map[v.idx] = new_edge
                self.tree.insert(new_edge, y)
                
            elif v.type == VertexType.REGULAR_RIGHT:
                left_edge = self.tree.find_left_edge(v.x, y)
                if left_edge:
                    if self.V[left_edge.helper_idx].type == VertexType.MERGE:
                        self.diagonals.append((v.idx, left_edge.helper_idx))
                    left_edge.helper_idx = v.idx

        t2 = time.time()
        self.stats['sweep_time'] = t2 - t1
        self.stats['splay_ops'] = self.tree.splay_ops
        self.stats['comparisons'] = self.tree.comparisons
        return self.diagonals, self.stats

# ============================================================
# Benchmark
# ============================================================

def create_comb_jitter(n, seed=42):
    random.seed(seed)
    teeth = n // 4
    pts = []
    # Start bottom left
    pts.append((-10, -50))
    for i in range(teeth):
        x_base = i * 4
        pts.append((x_base, random.random() * 20)) 
        pts.append((x_base + 1, 150 + random.random())) 
        pts.append((x_base + 2, random.random() * 20)) 
        pts.append((x_base + 3, random.random() * 10)) 
    pts.append((teeth*4 + 10, -50))
    return pts

def create_convex(n):
    return [(math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)) for i in range(n)]

def benchmark():
    print("=" * 80)
    print("LINEAR TRIANGULATION (BUCKET SORT + SPLAY) BENCHMARK")
    print("=" * 80)
    print(f"{'Type':<10} {'n':>6} {'Rotations':>10} {'Rot/N':>8} {'Comps':>10} {'Cmp/N':>8} {'Total(s)':>10}")
    
    cases = [
        ("Jitter", create_comb_jitter),
        ("Convex", create_convex)
    ]

    for name, gen in cases:
        print(f"--- {name} ---")
        for n in [1000, 5000, 10000, 20000, 50000]:
            pts = gen(n)
            tri = LinearTriangulator(pts)
            _, stats = tri.triangulate()
            
            ops = stats['splay_ops']
            cmps = stats['comparisons']
            total_t = stats['sort_time'] + stats['sweep_time']
            
            print(f"{name:<10} {n:>6} {ops:>10} {ops/n:>8.2f} {cmps:>10} {cmps/n:>8.2f} {total_t:>10.4f}")

if __name__ == "__main__":
    benchmark()

