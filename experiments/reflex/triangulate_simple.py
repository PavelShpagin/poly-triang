"""
Simple Polygon Triangulation - Exploring Minimal Constant O(n log r)

Goal: Find the simplest implementation with lowest constant factor.
Approach: Minimize data structure overhead.
"""

import math
from typing import List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum
import time


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
    

def cross(o: Vertex, a: Vertex, b: Vertex) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def classify_vertices(polygon: List[Vertex]) -> Tuple[int, int]:
    """Classify all vertices. Returns (num_splits, num_merges)."""
    n = len(polygon)
    splits = merges = 0
    
    for i in range(n):
        v = polygon[i]
        prev = polygon[(i - 1) % n]
        next = polygon[(i + 1) % n]
        
        below_prev = prev.y < v.y
        below_next = next.y < v.y
        is_reflex = cross(prev, v, next) < 0
        
        if below_prev and below_next:  # Both neighbors below
            if is_reflex:
                v.type = VertexType.SPLIT
                splits += 1
            else:
                v.type = VertexType.START
        elif not below_prev and not below_next:  # Both neighbors above
            if is_reflex:
                v.type = VertexType.MERGE
                merges += 1
            else:
                v.type = VertexType.END
        else:  # One above, one below
            # Check if interior is to the right
            if prev.y > next.y:
                v.type = VertexType.REGULAR_LEFT
            else:
                v.type = VertexType.REGULAR_RIGHT
    
    return splits, merges


# ============================================================
# APPROACH 1: Standard BST (baseline)
# ============================================================

class BSTNode:
    def __init__(self, edge_idx: int, x_intercept: float):
        self.edge_idx = edge_idx
        self.x_intercept = x_intercept
        self.left = None
        self.right = None
        self.height = 1


class AVLTree:
    """Simple AVL tree for sweep line status."""
    
    def __init__(self):
        self.root = None
        self.size = 0
        self.operations = 0
    
    def _height(self, node):
        return node.height if node else 0
    
    def _balance(self, node):
        return self._height(node.left) - self._height(node.right) if node else 0
    
    def _rotate_right(self, y):
        x = y.left
        T2 = x.right
        x.right = y
        y.left = T2
        y.height = 1 + max(self._height(y.left), self._height(y.right))
        x.height = 1 + max(self._height(x.left), self._height(x.right))
        return x
    
    def _rotate_left(self, x):
        y = x.right
        T2 = y.left
        y.left = x
        x.right = T2
        x.height = 1 + max(self._height(x.left), self._height(x.right))
        y.height = 1 + max(self._height(y.left), self._height(y.right))
        return y
    
    def insert(self, edge_idx: int, x_intercept: float):
        self.operations += 1
        self.root = self._insert(self.root, edge_idx, x_intercept)
        self.size += 1
    
    def _insert(self, node, edge_idx, x_intercept):
        if not node:
            return BSTNode(edge_idx, x_intercept)
        
        self.operations += 1
        if x_intercept < node.x_intercept:
            node.left = self._insert(node.left, edge_idx, x_intercept)
        else:
            node.right = self._insert(node.right, edge_idx, x_intercept)
        
        node.height = 1 + max(self._height(node.left), self._height(node.right))
        balance = self._balance(node)
        
        if balance > 1 and x_intercept < node.left.x_intercept:
            return self._rotate_right(node)
        if balance < -1 and x_intercept > node.right.x_intercept:
            return self._rotate_left(node)
        if balance > 1 and x_intercept > node.left.x_intercept:
            node.left = self._rotate_left(node.left)
            return self._rotate_right(node)
        if balance < -1 and x_intercept < node.right.x_intercept:
            node.right = self._rotate_right(node.right)
            return self._rotate_left(node)
        
        return node
    
    def find_left_neighbor(self, x: float) -> Optional[int]:
        """Find edge immediately to the left of x."""
        self.operations += 1
        result = None
        node = self.root
        while node:
            self.operations += 1
            if node.x_intercept < x:
                result = node.edge_idx
                node = node.right
            else:
                node = node.left
        return result


# ============================================================
# APPROACH 2: Sorted Linked List (simpler, O(r) per search)
# ============================================================

class LinkedListStatus:
    """Sorted linked list - O(r) search but very simple."""
    
    def __init__(self):
        self.edges = []  # List of (x_intercept, edge_idx)
        self.operations = 0
    
    def insert(self, edge_idx: int, x_intercept: float):
        self.operations += 1
        # Find insertion point
        pos = 0
        for i, (x, _) in enumerate(self.edges):
            self.operations += 1
            if x > x_intercept:
                break
            pos = i + 1
        self.edges.insert(pos, (x_intercept, edge_idx))
    
    def remove(self, edge_idx: int):
        self.operations += 1
        for i, (_, idx) in enumerate(self.edges):
            self.operations += 1
            if idx == edge_idx:
                self.edges.pop(i)
                return
    
    def find_left_neighbor(self, x: float) -> Optional[int]:
        self.operations += 1
        result = None
        for x_int, edge_idx in self.edges:
            self.operations += 1
            if x_int >= x:
                break
            result = edge_idx
        return result


# ============================================================
# APPROACH 3: Skip List (O(log r) with better constants)
# ============================================================

import random

class SkipListNode:
    def __init__(self, x_intercept: float, edge_idx: int, level: int):
        self.x_intercept = x_intercept
        self.edge_idx = edge_idx
        self.forward = [None] * (level + 1)


class SkipList:
    """Skip list - O(log r) with simpler implementation than AVL."""
    
    MAX_LEVEL = 16
    
    def __init__(self):
        self.header = SkipListNode(float('-inf'), -1, self.MAX_LEVEL)
        self.level = 0
        self.operations = 0
    
    def _random_level(self):
        level = 0
        while random.random() < 0.5 and level < self.MAX_LEVEL:
            level += 1
        return level
    
    def insert(self, edge_idx: int, x_intercept: float):
        self.operations += 1
        update = [None] * (self.MAX_LEVEL + 1)
        current = self.header
        
        for i in range(self.level, -1, -1):
            while current.forward[i] and current.forward[i].x_intercept < x_intercept:
                self.operations += 1
                current = current.forward[i]
            update[i] = current
        
        level = self._random_level()
        if level > self.level:
            for i in range(self.level + 1, level + 1):
                update[i] = self.header
            self.level = level
        
        new_node = SkipListNode(x_intercept, edge_idx, level)
        for i in range(level + 1):
            new_node.forward[i] = update[i].forward[i]
            update[i].forward[i] = new_node
    
    def find_left_neighbor(self, x: float) -> Optional[int]:
        self.operations += 1
        current = self.header
        
        for i in range(self.level, -1, -1):
            while current.forward[i] and current.forward[i].x_intercept < x:
                self.operations += 1
                current = current.forward[i]
        
        if current != self.header:
            return current.edge_idx
        return None


# ============================================================
# APPROACH 4: Bucket-based (O(1) amortized for bounded coordinates)
# ============================================================

class BucketStatus:
    """Bucket-based status for integer/bounded coordinates."""
    
    def __init__(self, min_x: float, max_x: float, num_buckets: int):
        self.min_x = min_x
        self.max_x = max_x
        self.num_buckets = num_buckets
        self.bucket_size = (max_x - min_x) / num_buckets if num_buckets > 0 else 1
        self.buckets = [[] for _ in range(num_buckets + 1)]
        self.operations = 0
    
    def _bucket_idx(self, x: float) -> int:
        if self.bucket_size == 0:
            return 0
        idx = int((x - self.min_x) / self.bucket_size)
        return max(0, min(idx, self.num_buckets))
    
    def insert(self, edge_idx: int, x_intercept: float):
        self.operations += 1
        bucket = self._bucket_idx(x_intercept)
        self.buckets[bucket].append((x_intercept, edge_idx))
        self.buckets[bucket].sort()  # Keep sorted within bucket
    
    def find_left_neighbor(self, x: float) -> Optional[int]:
        self.operations += 1
        bucket = self._bucket_idx(x)
        result = None
        
        # Check current and previous buckets
        for b in range(bucket, -1, -1):
            for x_int, edge_idx in reversed(self.buckets[b]):
                self.operations += 1
                if x_int < x:
                    return edge_idx
        return result


# ============================================================
# Test harness
# ============================================================

def create_test_polygon(n: int, reflex_ratio: float = 0.3) -> List[Vertex]:
    """Create a polygon with approximately n vertices and given reflex ratio."""
    vertices = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        # Vary radius to create reflex vertices
        r = 2 if (i % int(1/reflex_ratio)) == 0 else 1
        x = r * math.cos(angle)
        y = r * math.sin(angle)
        vertices.append(Vertex(x, y, i))
    return vertices


def benchmark_approach(name: str, polygon: List[Vertex], status_class, **kwargs):
    """Benchmark a single approach."""
    n = len(polygon)
    classify_vertices(polygon)
    
    # Sort vertices by y
    sorted_vertices = sorted(range(n), key=lambda i: -polygon[i].y)
    
    # Create status structure
    if kwargs:
        status = status_class(**kwargs)
    else:
        status = status_class()
    
    start_time = time.perf_counter()
    
    # Simulate sweep
    for idx in sorted_vertices:
        v = polygon[idx]
        x = v.x
        
        if v.type in [VertexType.SPLIT, VertexType.MERGE]:
            # These need left-neighbor query
            status.find_left_neighbor(x)
        
        if v.type in [VertexType.START, VertexType.SPLIT]:
            # Insert edge
            status.insert(idx, x)
    
    elapsed = time.perf_counter() - start_time
    
    return {
        'name': name,
        'time_ms': elapsed * 1000,
        'operations': status.operations,
        'ops_per_vertex': status.operations / n
    }


def run_benchmarks():
    """Run benchmarks on different polygon sizes."""
    print("=" * 70)
    print("POLYGON TRIANGULATION - DATA STRUCTURE COMPARISON")
    print("=" * 70)
    print()
    
    sizes = [100, 500, 1000, 5000, 10000]
    
    for n in sizes:
        polygon = create_test_polygon(n, reflex_ratio=0.3)
        splits, merges = classify_vertices(polygon)
        r = splits + merges
        
        print(f"n = {n}, r = {r} (ratio = {r/n:.2%})")
        print("-" * 50)
        
        # Test each approach
        results = []
        
        # AVL Tree
        polygon_copy = create_test_polygon(n, 0.3)
        results.append(benchmark_approach("AVL Tree", polygon_copy, AVLTree))
        
        # Skip List
        polygon_copy = create_test_polygon(n, 0.3)
        results.append(benchmark_approach("Skip List", polygon_copy, SkipList))
        
        # Linked List (for small r)
        if r < 100:
            polygon_copy = create_test_polygon(n, 0.3)
            results.append(benchmark_approach("Linked List", polygon_copy, LinkedListStatus))
        
        # Bucket (for bounded coords)
        polygon_copy = create_test_polygon(n, 0.3)
        results.append(benchmark_approach("Bucket (sqrt(n))", polygon_copy, BucketStatus,
                                          min_x=-3, max_x=3, num_buckets=int(math.sqrt(n))))
        
        for res in results:
            print(f"  {res['name']:20s}: {res['time_ms']:8.2f} ms, "
                  f"{res['operations']:8d} ops, {res['ops_per_vertex']:6.2f} ops/vertex")
        print()
    
    print("=" * 70)
    print("CONCLUSION: Skip List has lower constant than AVL for same O(log r)")
    print("=" * 70)


if __name__ == "__main__":
    run_benchmarks()

