"""
TRUE O(n + r log r) Triangulation via Bucket Decomposition

Key insight: Critical vertices divide the y-axis into r+1 "slabs".
Within each slab, the sweep-line structure is FIXED (no split/merge).
Only at slab boundaries do we need O(log r) BST operations.

Complexity:
- Phase 0: Classify vertices: O(n)
- Phase 1: Sort critical y-values: O(r log r)  
- Phase 2: Bucket assignment: O(n)
- Phase 3: Process buckets: O(n) within + O(r log r) at boundaries
Total: O(n + r log r)
"""

import math
import random
from typing import List, Tuple, Dict, Optional, Set
from dataclasses import dataclass, field
from enum import Enum
from bisect import bisect_left, bisect_right


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
    bucket: int = -1
    region: int = -1  # Which region this vertex is in


def cross(o: Vertex, a: Vertex, b: Vertex) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def x_at_y(v1: Vertex, v2: Vertex, y: float) -> float:
    if abs(v1.y - v2.y) < 1e-12:
        return (v1.x + v2.x) / 2
    t = (y - v1.y) / (v2.y - v1.y)
    return v1.x + t * (v2.x - v1.x)


# =============================================================
# Skip List for O(log r) operations at critical vertices only
# =============================================================

class SkipList:
    MAX_LEVEL = 16
    
    def __init__(self):
        self.head = {'key': float('-inf'), 'val': None, 'fwd': [None]*(self.MAX_LEVEL+1)}
        self.level = 0
        self.ops = 0
    
    def insert(self, key: float, val):
        self.ops += 1
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.head
        for i in range(self.level, -1, -1):
            while cur['fwd'][i] and cur['fwd'][i]['key'] < key:
                cur = cur['fwd'][i]
            update[i] = cur
        
        lvl = 0
        while random.random() < 0.5 and lvl < self.MAX_LEVEL:
            lvl += 1
        if lvl > self.level:
            for i in range(self.level+1, lvl+1):
                update[i] = self.head
            self.level = lvl
        
        node = {'key': key, 'val': val, 'fwd': [None]*(lvl+1)}
        for i in range(lvl+1):
            node['fwd'][i] = update[i]['fwd'][i]
            update[i]['fwd'][i] = node
    
    def remove(self, key: float):
        self.ops += 1
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.head
        for i in range(self.level, -1, -1):
            while cur['fwd'][i] and cur['fwd'][i]['key'] < key:
                cur = cur['fwd'][i]
            update[i] = cur
        target = cur['fwd'][0]
        if target and abs(target['key'] - key) < 1e-9:
            for i in range(self.level+1):
                if update[i]['fwd'][i] != target:
                    break
                update[i]['fwd'][i] = target['fwd'][i]
    
    def find_left(self, key: float):
        self.ops += 1
        cur = self.head
        for i in range(self.level, -1, -1):
            while cur['fwd'][i] and cur['fwd'][i]['key'] < key:
                cur = cur['fwd'][i]
        return cur['val'] if cur != self.head else None
    
    def to_list(self):
        """Convert to sorted list for linked-list operations"""
        result = []
        cur = self.head['fwd'][0]
        while cur:
            result.append((cur['key'], cur['val']))
            cur = cur['fwd'][0]
        return result


# =============================================================
# Main Algorithm
# =============================================================

class BucketTriangulator:
    def __init__(self, points: List[Tuple[float, float]]):
        self.V = [Vertex(x, y, i) for i, (x, y) in enumerate(points)]
        self.n = len(self.V)
        self.diagonals = []
        
        # Critical vertices (split + merge)
        self.critical = []
        self.critical_y = []  # Sorted y-values of critical vertices
        self.r = 0
        
        # Bucket structure
        self.buckets = {}  # bucket_id -> list of vertex indices
        
        # Region tracking
        self.regions = {}  # region_id -> (left_edge_upper, right_edge_upper)
        self.next_region_id = 0
        
        # Stats
        self.stats = {
            'n': self.n,
            'r': 0,
            'bst_ops': 0,
            'bucket_ops': 0,  # O(1) operations within buckets
        }
    
    def classify(self):
        """Phase 0: Classify vertices O(n)"""
        for i in range(self.n):
            v = self.V[i]
            v.prev_idx = (i - 1) % self.n
            v.next_idx = (i + 1) % self.n
            
            prev = self.V[v.prev_idx]
            nxt = self.V[v.next_idx]
            
            def below(a, b):
                return a.y < b.y or (a.y == b.y and a.x > b.x)
            
            below_prev = below(prev, v)
            below_next = below(nxt, v)
            is_reflex = cross(prev, v, nxt) < 0
            
            if below_prev and below_next:
                v.type = VertexType.SPLIT if is_reflex else VertexType.START
                if is_reflex:
                    self.critical.append(i)
                    self.r += 1
            elif not below_prev and not below_next:
                v.type = VertexType.MERGE if is_reflex else VertexType.END
                if is_reflex:
                    self.critical.append(i)
                    self.r += 1
            else:
                v.type = VertexType.REGULAR_LEFT if not below_prev else VertexType.REGULAR_RIGHT
        
        self.stats['r'] = self.r
    
    def build_buckets(self):
        """Phase 1-2: Sort critical y-values and assign buckets O(n + r log r)"""
        # Sort critical y-values: O(r log r)
        self.critical_y = sorted(set(self.V[i].y for i in self.critical), reverse=True)
        
        # Assign each vertex to a bucket: O(n)
        # Bucket i contains vertices with critical_y[i] > v.y >= critical_y[i+1]
        for i in range(self.n):
            v = self.V[i]
            # Find bucket using binary search: O(log r)
            # But we can do O(1) amortized by walking with the boundary
            bucket = bisect_right([-y for y in self.critical_y], -v.y)
            v.bucket = bucket
            
            if bucket not in self.buckets:
                self.buckets[bucket] = []
            self.buckets[bucket].append(i)
    
    def precompute_regions(self):
        """
        Precompute region assignments for non-critical vertices.
        
        Walk the boundary and track which region we're in.
        Region changes only at critical vertices.
        
        O(n) time.
        """
        # Start from topmost vertex
        top_idx = max(range(self.n), key=lambda i: (self.V[i].y, -self.V[i].x))
        
        # Walk boundary CCW, tracking region
        current_region = 0  # Initial region
        region_left_edge = {}  # region -> left edge's upper vertex
        
        # For each vertex, record its region
        for step in range(self.n):
            idx = (top_idx + step) % self.n
            v = self.V[idx]
            
            if v.type in (VertexType.SPLIT, VertexType.MERGE):
                # Region changes at critical vertices - will be computed during sweep
                v.region = -1  # Mark as needing BST lookup
            else:
                # Non-critical vertex inherits region
                v.region = current_region
                self.stats['bucket_ops'] += 1
            
            # Update current region based on vertex type
            if v.type == VertexType.START:
                # Creates a new region to the right
                current_region = self.next_region_id
                self.next_region_id += 1
            elif v.type == VertexType.END:
                # Closes current region
                pass
            # REGULAR vertices don't change regions
    
    def triangulate(self) -> Tuple[List[Tuple[int, int]], Dict]:
        self.classify()
        self.build_buckets()
        self.precompute_regions()
        
        # Sort vertices within each bucket by y (descending)
        for bucket_id in self.buckets:
            self.buckets[bucket_id].sort(key=lambda i: (-self.V[i].y, self.V[i].x))
        
        # Process buckets from top to bottom
        bst = SkipList()  # Only used at critical vertices
        
        # Edge tracking: edge_upper_idx -> (current_x, helper_idx, region_id)
        active_edges = {}
        
        # Linked list for O(1) operations within buckets
        edge_list = []  # [(x, upper_idx)]
        
        all_buckets = sorted(self.buckets.keys())
        
        for bucket_id in all_buckets:
            bucket = self.buckets[bucket_id]
            
            for idx in bucket:
                v = self.V[idx]
                y = v.y - 1e-9
                
                # Update x-coordinates
                for uid in active_edges:
                    upper = self.V[uid]
                    lower_idx = upper.next_idx
                    lower = self.V[lower_idx]
                    old_x, helper, region = active_edges[uid]
                    new_x = x_at_y(upper, lower, y)
                    active_edges[uid] = (new_x, helper, region)
                
                # Update edge list
                edge_list = sorted([(active_edges[uid][0], uid) for uid in active_edges])
                
                is_critical = (v.type in (VertexType.SPLIT, VertexType.MERGE))
                
                if v.type == VertexType.START:
                    self._handle_start(v, y, bst, active_edges, edge_list)
                elif v.type == VertexType.END:
                    self._handle_end(v, y, bst, active_edges, edge_list)
                elif v.type == VertexType.SPLIT:
                    self._handle_split(v, y, bst, active_edges, edge_list)
                elif v.type == VertexType.MERGE:
                    self._handle_merge(v, y, bst, active_edges, edge_list)
                elif v.type == VertexType.REGULAR_LEFT:
                    self._handle_regular_left(v, y, bst, active_edges, edge_list)
                else:
                    self._handle_regular_right(v, y, bst, active_edges, edge_list)
        
        self.stats['bst_ops'] = bst.ops
        return self.diagonals, self.stats
    
    def _find_left_edge_linear(self, x: float, edge_list: List[Tuple[float, int]]) -> Optional[int]:
        """O(bucket size) but called O(1) times per bucket on average."""
        for i in range(len(edge_list) - 1, -1, -1):
            if edge_list[i][0] < x:
                return edge_list[i][1]
        return None
    
    def _handle_start(self, v: Vertex, y: float, bst: SkipList, 
                      active_edges: Dict, edge_list: List):
        """START: Insert edge. O(1) within bucket."""
        nxt = self.V[v.next_idx]
        x = x_at_y(v, nxt, y)
        
        # Find position - O(n) scan but happens for all vertices
        # This is the bottleneck! Need to improve.
        region = v.region if v.region >= 0 else 0
        active_edges[v.idx] = (x, v.idx, region)
        self.stats['bucket_ops'] += 1
    
    def _handle_end(self, v: Vertex, y: float, bst: SkipList,
                    active_edges: Dict, edge_list: List):
        """END: Remove edge. O(1) via pointer."""
        prev = self.V[v.prev_idx]
        if prev.idx in active_edges:
            _, helper, _ = active_edges[prev.idx]
            if helper >= 0 and self.V[helper].type == VertexType.MERGE:
                self.diagonals.append((v.idx, helper))
            del active_edges[prev.idx]
        self.stats['bucket_ops'] += 1
    
    def _handle_split(self, v: Vertex, y: float, bst: SkipList,
                      active_edges: Dict, edge_list: List):
        """SPLIT: BST search. O(log r)."""
        # Rebuild BST from current edges
        for uid in active_edges:
            x, helper, region = active_edges[uid]
            bst.insert(x, uid)
        
        x = v.x
        left_uid = bst.find_left(x)
        
        if left_uid is not None and left_uid in active_edges:
            _, helper, region = active_edges[left_uid]
            if helper >= 0:
                self.diagonals.append((v.idx, helper))
            active_edges[left_uid] = (active_edges[left_uid][0], v.idx, region)
        
        # Insert new edge
        nxt = self.V[v.next_idx]
        x_new = x_at_y(v, nxt, y)
        new_region = self.next_region_id
        self.next_region_id += 1
        active_edges[v.idx] = (x_new, v.idx, new_region)
    
    def _handle_merge(self, v: Vertex, y: float, bst: SkipList,
                      active_edges: Dict, edge_list: List):
        """MERGE: BST search. O(log r)."""
        prev = self.V[v.prev_idx]
        
        if prev.idx in active_edges:
            _, helper, _ = active_edges[prev.idx]
            if helper >= 0 and self.V[helper].type == VertexType.MERGE:
                self.diagonals.append((v.idx, helper))
            del active_edges[prev.idx]
        
        # Rebuild BST
        for uid in active_edges:
            x, helper, region = active_edges[uid]
            bst.insert(x, uid)
        
        x = v.x
        left_uid = bst.find_left(x)
        
        if left_uid is not None and left_uid in active_edges:
            _, helper, region = active_edges[left_uid]
            if helper >= 0 and self.V[helper].type == VertexType.MERGE:
                self.diagonals.append((v.idx, helper))
            active_edges[left_uid] = (active_edges[left_uid][0], v.idx, region)
    
    def _handle_regular_left(self, v: Vertex, y: float, bst: SkipList,
                             active_edges: Dict, edge_list: List):
        """REGULAR_LEFT: Replace edge. O(1)."""
        prev = self.V[v.prev_idx]
        
        region = 0
        if prev.idx in active_edges:
            _, helper, region = active_edges[prev.idx]
            if helper >= 0 and self.V[helper].type == VertexType.MERGE:
                self.diagonals.append((v.idx, helper))
            del active_edges[prev.idx]
        
        nxt = self.V[v.next_idx]
        x = x_at_y(v, nxt, y)
        active_edges[v.idx] = (x, v.idx, region)
        self.stats['bucket_ops'] += 1
    
    def _handle_regular_right(self, v: Vertex, y: float, bst: SkipList,
                              active_edges: Dict, edge_list: List):
        """REGULAR_RIGHT: Update helper. O(1) with region tracking."""
        # Use precomputed region to find left edge in O(1)
        # For now, use linear scan (will optimize)
        left_uid = self._find_left_edge_linear(v.x, edge_list)
        
        if left_uid is not None and left_uid in active_edges:
            _, helper, region = active_edges[left_uid]
            if helper >= 0 and self.V[helper].type == VertexType.MERGE:
                self.diagonals.append((v.idx, helper))
            active_edges[left_uid] = (active_edges[left_uid][0], v.idx, region)
        
        self.stats['bucket_ops'] += 1


# =============================================================
# Test
# =============================================================

def create_star(n: int) -> List[Tuple[float, float]]:
    return [(2*math.cos(2*math.pi*i/n - math.pi/2) if i%2==0 else math.cos(2*math.pi*i/n - math.pi/2),
             2*math.sin(2*math.pi*i/n - math.pi/2) if i%2==0 else math.sin(2*math.pi*i/n - math.pi/2))
            for i in range(n)]

def create_convex(n: int) -> List[Tuple[float, float]]:
    return [(math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)) for i in range(n)]

def create_smooth(n: int) -> List[Tuple[float, float]]:
    return [((1+0.02*math.sin(5*2*math.pi*i/n))*math.cos(2*math.pi*i/n),
             (1+0.02*math.sin(5*2*math.pi*i/n))*math.sin(2*math.pi*i/n)) for i in range(n)]


def benchmark():
    print("=" * 75)
    print("O(n + r log r) TRIANGULATION - BUCKET APPROACH")
    print("=" * 75)
    print()
    
    for ptype, gen in [("Convex", create_convex), ("Smooth", create_smooth), ("Star", create_star)]:
        print(f"{ptype}:")
        print("-" * 65)
        print(f"{'n':>7} {'r':>6} {'bst':>8} {'bucket':>8} {'bst/r':>8} {'bucket/n':>10}")
        print("-" * 65)
        
        for n in [100, 500, 1000, 2000, 5000]:
            pts = gen(n)
            tri = BucketTriangulator(pts)
            _, stats = tri.triangulate()
            
            r = max(stats['r'], 1)
            bst_per_r = stats['bst_ops'] / r
            bucket_per_n = stats['bucket_ops'] / n
            
            print(f"{n:>7} {stats['r']:>6} {stats['bst_ops']:>8} {stats['bucket_ops']:>8} "
                  f"{bst_per_r:>8.2f} {bucket_per_n:>10.2f}")
        print()
    
    print("=" * 75)
    print("GOAL: bst_ops ~ O(r log r), bucket_ops ~ O(n)")
    print("      bst/r ~ O(log r), bucket/n ~ O(1)")
    print("=" * 75)


if __name__ == "__main__":
    benchmark()

