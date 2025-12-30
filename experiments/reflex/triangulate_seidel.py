
import random
import sys
from dataclasses import dataclass
from typing import List, Optional

sys.setrecursionlimit(20000)

@dataclass
class Point:
    x: float
    y: float
    id: int

    def __hash__(self):
        return self.id
    def __eq__(self, other):
        if not isinstance(other, Point): return False
        return self.id == other.id

@dataclass
class Segment:
    p: Point
    q: Point
    id: int
    
    def __post_init__(self):
        if (self.p.x > self.q.x) or (self.p.x == self.q.x and self.p.y > self.q.y):
            self.p, self.q = self.q, self.p

class Trapezoid:
    def __init__(self, top, bottom, left_p, right_p):
        self.top = top        # Segment
        self.bottom = bottom  # Segment
        self.left_p = left_p  # Point
        self.right_p = right_p# Point
        
        # Neighbors
        self.u_left = None
        self.l_left = None
        self.u_right = None
        self.l_right = None
        
        # DAG Node
        self.node = None
        
        self.is_inside = False

class Node:
    def __init__(self, type_str, item, left=None, right=None):
        self.type = type_str # 'x', 'y', 'leaf'
        self.item = item
        self.left = left
        self.right = right
        self.parents = []

def is_left_of_segment(seg, p):
    # Cross product
    return (seg.q.x - seg.p.x) * (p.y - seg.p.y) - (seg.q.y - seg.p.y) * (p.x - seg.p.x) > 0

class SeidelTriangulator:
    def __init__(self, points: List[tuple]):
        self.points = [Point(x, y, i) for i, (x, y) in enumerate(points)]
        self.n = len(points)
        self.segments = []
        for i in range(self.n):
            self.segments.append(Segment(self.points[i], self.points[(i+1)%self.n], i))
        
        # Shuffle segments for randomization
        self.segments_shuffled = self.segments[:]
        random.shuffle(self.segments_shuffled)
        
        # Init structure
        min_x = min(p.x for p in self.points) - 1
        max_x = max(p.x for p in self.points) + 1
        min_y = min(p.y for p in self.points) - 1
        max_y = max(p.y for p in self.points) + 1
        
        self.p_min = Point(min_x, min_y, -1)
        self.p_max = Point(max_x, max_y, -2)
        
        top = Segment(Point(min_x, max_y, -3), Point(max_x, max_y, -4), -1)
        bottom = Segment(Point(min_x, min_y, -5), Point(max_x, min_y, -6), -2)
        
        self.root_trap = Trapezoid(top, bottom, self.p_min, self.p_max)
        self.dag_root = Node('leaf', self.root_trap)
        self.root_trap.node = self.dag_root
        
        self.trapezoids = {self.root_trap}
        self.ops = 0

    def query(self, p, start_node=None):
        node = start_node if start_node else self.dag_root
        while node.type != 'leaf':
            self.ops += 1
            if node.type == 'x':
                if p.x == node.item.x:
                    # Point equality/coincidence handling is tricky.
                    # Assume p is slightly right.
                    node = node.right
                elif p.x < node.item.x:
                    node = node.left
                else:
                    node = node.right
            elif node.type == 'y':
                if is_left_of_segment(node.item, p):
                    node = node.left
                else:
                    node = node.right
        return node.item

    def update_dag(self, old_trap, new_subtree):
        # Replace old_trap's leaf node with new_subtree in parents
        leaf_node = old_trap.node
        for parent in leaf_node.parents:
            if parent.left == leaf_node:
                parent.left = new_subtree
            else:
                parent.right = new_subtree
        
        # If root was replaced
        if leaf_node == self.dag_root:
            self.dag_root = new_subtree
            
    def insert_segment(self, seg):
        # Locate left endpoint
        t_first = self.query(seg.p)
        t_last = self.query(seg.q)
        
        # Traverse trapezoids
        # This is the hard part: finding ALL intersected trapezoids
        # Requires following neighbor pointers
        pass

# ...
# Implementing Seidel is too complex for a single-file prototype in < 5 mins.
# Let's pivot to "Ear Clipping with Geometric Hashing" (FIST).
# It's robust, simpler, and "Industrial Strength".
# It might not be O(n) worst case, but O(n) practical.
# Wait, user said "reach holy grail".
# Seidel IS the holy grail.

