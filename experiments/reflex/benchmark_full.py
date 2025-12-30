"""
Full benchmark for polygon triangulation data structures.
Tests various polygon types and sizes.
"""

import math
import random
import time
from dataclasses import dataclass
from typing import List, Tuple
from enum import Enum


class VertexType(Enum):
    START = 1
    END = 2
    SPLIT = 3
    MERGE = 4
    REGULAR = 5


@dataclass
class Vertex:
    x: float
    y: float
    idx: int
    type: VertexType = None


def cross(o: Vertex, a: Vertex, b: Vertex) -> float:
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def classify(polygon: List[Vertex]) -> int:
    """Classify vertices, return r."""
    n = len(polygon)
    r = 0
    for i in range(n):
        v = polygon[i]
        prev = polygon[(i - 1) % n]
        nxt = polygon[(i + 1) % n]
        
        below_prev = prev.y < v.y
        below_next = nxt.y < v.y
        is_reflex = cross(prev, v, nxt) < 0
        
        if below_prev and below_next:
            v.type = VertexType.SPLIT if is_reflex else VertexType.START
            if is_reflex: r += 1
        elif not below_prev and not below_next:
            v.type = VertexType.MERGE if is_reflex else VertexType.END
            if is_reflex: r += 1
        else:
            v.type = VertexType.REGULAR
    return r


# ============================================================
# Simple Skip List Implementation
# ============================================================

class SkipNode:
    def __init__(self, key: float, val: int, level: int):
        self.key = key
        self.val = val
        self.forward = [None] * (level + 1)


class SkipList:
    def __init__(self, max_level: int = 16):
        self.max_level = max_level
        self.header = SkipNode(float('-inf'), -1, max_level)
        self.level = 0
        self.ops = 0
    
    def _rand_level(self) -> int:
        lvl = 0
        while random.random() < 0.5 and lvl < self.max_level:
            lvl += 1
        return lvl
    
    def insert(self, key: float, val: int):
        update = [None] * (self.max_level + 1)
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
        
        node = SkipNode(key, val, lvl)
        for i in range(lvl + 1):
            node.forward[i] = update[i].forward[i]
            update[i].forward[i] = node
    
    def delete(self, key: float):
        update = [None] * (self.max_level + 1)
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
    
    def find_left(self, key: float) -> int:
        cur = self.header
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                self.ops += 1
                cur = cur.forward[i]
        return cur.val if cur != self.header else -1


# ============================================================
# Test Polygon Generators
# ============================================================

def star_polygon(n: int) -> List[Vertex]:
    """Star with alternating radii - high reflex ratio."""
    verts = []
    for i in range(n):
        angle = 2 * math.pi * i / n - math.pi/2
        r = 2 if i % 2 == 0 else 1
        verts.append(Vertex(r * math.cos(angle), r * math.sin(angle), i))
    return verts


def smooth_polygon(n: int) -> List[Vertex]:
    """Nearly convex - low reflex ratio."""
    verts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 1 + 0.1 * math.sin(5 * angle)  # Small perturbation
        verts.append(Vertex(r * math.cos(angle), r * math.sin(angle), i))
    return verts


def convex_polygon(n: int) -> List[Vertex]:
    """Convex - r = 0."""
    verts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        verts.append(Vertex(math.cos(angle), math.sin(angle), i))
    return verts


def zigzag_polygon(n: int) -> List[Vertex]:
    """Zigzag - many reflex vertices."""
    verts = []
    for i in range(n):
        x = i
        y = 0 if i % 2 == 0 else 1
        verts.append(Vertex(x, y, i))
    # Close with top edge
    return verts


# ============================================================
# Sweep Simulation
# ============================================================

def simulate_sweep(polygon: List[Vertex]) -> Tuple[int, float]:
    """Simulate sweep, return (operations, time)."""
    n = len(polygon)
    r = classify(polygon)
    
    # Sort by y
    order = sorted(range(n), key=lambda i: -polygon[i].y)
    
    status = SkipList()
    helper = {}  # edge_idx -> helper vertex
    diagonals = []
    
    start = time.perf_counter()
    
    for idx in order:
        v = polygon[idx]
        prev_idx = (idx - 1) % n
        next_idx = (idx + 1) % n
        
        if v.type == VertexType.START:
            edge_key = v.x
            status.insert(edge_key, idx)
            helper[idx] = idx
            
        elif v.type == VertexType.END:
            edge_key = polygon[prev_idx].x
            if prev_idx in helper and polygon[helper[prev_idx]].type == VertexType.MERGE:
                diagonals.append((idx, helper[prev_idx]))
            status.delete(edge_key)
            
        elif v.type == VertexType.SPLIT:
            left_edge = status.find_left(v.x)
            if left_edge >= 0 and left_edge in helper:
                diagonals.append((idx, helper[left_edge]))
                helper[left_edge] = idx
            edge_key = v.x
            status.insert(edge_key, idx)
            helper[idx] = idx
            
        elif v.type == VertexType.MERGE:
            edge_key = polygon[prev_idx].x
            if prev_idx in helper and polygon[helper[prev_idx]].type == VertexType.MERGE:
                diagonals.append((idx, helper[prev_idx]))
            status.delete(edge_key)
            
            left_edge = status.find_left(v.x)
            if left_edge >= 0 and left_edge in helper:
                if polygon[helper[left_edge]].type == VertexType.MERGE:
                    diagonals.append((idx, helper[left_edge]))
                helper[left_edge] = idx
                
        else:  # REGULAR
            # Simplified handling
            pass
    
    elapsed = time.perf_counter() - start
    return status.ops, elapsed


# ============================================================
# Main Benchmark
# ============================================================

def main():
    print("=" * 75)
    print("POLYGON TRIANGULATION BENCHMARK")
    print("=" * 75)
    print()
    
    # Test different polygon types
    test_cases = [
        ("Convex", convex_polygon),
        ("Smooth", smooth_polygon),
        ("Star", star_polygon),
    ]
    
    sizes = [100, 500, 1000, 2000, 5000, 10000]
    
    for name, generator in test_cases:
        print(f"\n{name} Polygons:")
        print("-" * 60)
        print(f"{'n':>8} {'r':>8} {'r/n':>8} {'ops':>10} {'ops/n':>8} {'ms':>8}")
        print("-" * 60)
        
        for n in sizes:
            polygon = generator(n)
            r = classify(polygon)
            
            # Run multiple times for stability
            total_ops = 0
            total_time = 0
            runs = 3
            
            for _ in range(runs):
                polygon = generator(n)
                classify(polygon)
                ops, elapsed = simulate_sweep(polygon)
                total_ops += ops
                total_time += elapsed
            
            avg_ops = total_ops / runs
            avg_time = total_time / runs * 1000
            
            print(f"{n:>8} {r:>8} {r/n:>8.2%} {avg_ops:>10.0f} {avg_ops/n:>8.2f} {avg_time:>8.2f}")
    
    print()
    print("=" * 75)
    print("VERIFICATION: ops/n should grow as O(log r)")
    print("For Star polygons: r ~ n/2, so log(r) ~ log(n) - 1")
    print("For Smooth polygons: r ~ constant, so log(r) ~ constant")
    print("For Convex polygons: r = 0, so ops should be O(1)")
    print("=" * 75)


if __name__ == "__main__":
    main()

