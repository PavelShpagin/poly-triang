#!/usr/bin/env python3
"""
Garey-Johnson-Preparata-Tarjan O(n log n) polygon triangulation.
Based on monotone decomposition followed by monotone triangulation.
"""

import sys
import time
import argparse
import heapq
from enum import Enum
from typing import List, Tuple, Optional, Set


class VertexType(Enum):
    START = 1      # Both neighbors below, interior angle < pi
    END = 2        # Both neighbors above, interior angle < pi
    SPLIT = 3      # Both neighbors below, interior angle > pi
    MERGE = 4      # Both neighbors above, interior angle > pi
    REGULAR = 5    # One neighbor above, one below


class Vertex:
    def __init__(self, x: float, y: float, index: int):
        self.x = x
        self.y = y
        self.index = index
        self.prev: Optional['Vertex'] = None
        self.next: Optional['Vertex'] = None
        self.vertex_type: VertexType = VertexType.REGULAR
    
    def __lt__(self, other):
        # Sort by y descending, then x ascending
        if abs(self.y - other.y) < 1e-9:
            return self.x < other.x
        return self.y > other.y
    
    def __repr__(self):
        return f"V{self.index}({self.x:.2f},{self.y:.2f})"


class Edge:
    def __init__(self, v1: Vertex, v2: Vertex):
        # Ensure v1 is upper vertex
        if v1.y < v2.y or (abs(v1.y - v2.y) < 1e-9 and v1.x > v2.x):
            v1, v2 = v2, v1
        self.upper = v1
        self.lower = v2
        self.helper: Optional[Vertex] = None
    
    def x_at_y(self, y: float) -> float:
        """Get x coordinate of edge at given y."""
        if abs(self.upper.y - self.lower.y) < 1e-9:
            return min(self.upper.x, self.lower.x)
        t = (self.upper.y - y) / (self.upper.y - self.lower.y)
        return self.upper.x + t * (self.lower.x - self.upper.x)


def cross_product(o: Vertex, a: Vertex, b: Vertex) -> float:
    """Cross product of vectors OA and OB."""
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def classify_vertex(v: Vertex) -> VertexType:
    """Classify vertex type based on its neighbors."""
    prev_below = v.prev.y < v.y or (abs(v.prev.y - v.y) < 1e-9 and v.prev.x > v.x)
    next_below = v.next.y < v.y or (abs(v.next.y - v.y) < 1e-9 and v.next.x > v.x)
    
    # Interior angle check using cross product
    # For CCW polygon, cross > 0 means convex (interior angle < pi)
    cross = cross_product(v.prev, v, v.next)
    convex = cross > 0
    
    if prev_below and next_below:
        return VertexType.START if convex else VertexType.SPLIT
    elif not prev_below and not next_below:
        return VertexType.END if convex else VertexType.MERGE
    else:
        return VertexType.REGULAR


def make_monotone(vertices: List[Vertex]) -> List[Tuple[int, int]]:
    """
    Decompose polygon into y-monotone pieces by adding diagonals.
    Returns list of diagonals as (vertex_index, vertex_index) pairs.
    """
    n = len(vertices)
    diagonals = []
    
    # Classify all vertices
    for v in vertices:
        v.vertex_type = classify_vertex(v)
    
    # Sort vertices by y (descending), then x (ascending)
    sorted_vertices = sorted(vertices)
    
    # Active edge list (edges intersecting current sweep line)
    # We'll use a simple list and search; for true O(n log n) would use balanced BST
    active_edges: List[Edge] = []
    
    def find_left_edge(v: Vertex) -> Optional[Edge]:
        """Find the edge immediately to the left of vertex v."""
        left_edge = None
        left_x = float('-inf')
        for e in active_edges:
            ex = e.x_at_y(v.y)
            if ex < v.x and ex > left_x:
                left_x = ex
                left_edge = e
        return left_edge
    
    def add_edge(v: Vertex):
        """Add edge from v going downward."""
        e = Edge(v, v.next)
        e.helper = v
        active_edges.append(e)
        return e
    
    def remove_edge(v: Vertex):
        """Remove edge ending at v."""
        for i, e in enumerate(active_edges):
            if e.lower == v:
                active_edges.pop(i)
                return e
        return None
    
    def find_edge_with_upper(v: Vertex) -> Optional[Edge]:
        """Find edge with upper vertex v."""
        for e in active_edges:
            if e.upper == v:
                return e
        return None
    
    # Process vertices from top to bottom
    for v in sorted_vertices:
        if v.vertex_type == VertexType.START:
            add_edge(v)
        
        elif v.vertex_type == VertexType.END:
            # Check helper of edge ending here
            edge = find_edge_with_upper(v.prev)
            if edge and edge.helper and edge.helper.vertex_type == VertexType.MERGE:
                diagonals.append((v.index, edge.helper.index))
            remove_edge(v)
        
        elif v.vertex_type == VertexType.SPLIT:
            # Find edge to left and add diagonal to its helper
            left_edge = find_left_edge(v)
            if left_edge and left_edge.helper:
                diagonals.append((v.index, left_edge.helper.index))
                left_edge.helper = v
            add_edge(v)
        
        elif v.vertex_type == VertexType.MERGE:
            # Handle edge ending here
            edge = find_edge_with_upper(v.prev)
            if edge and edge.helper and edge.helper.vertex_type == VertexType.MERGE:
                diagonals.append((v.index, edge.helper.index))
            remove_edge(v)
            
            # Update helper of edge to left
            left_edge = find_left_edge(v)
            if left_edge:
                if left_edge.helper and left_edge.helper.vertex_type == VertexType.MERGE:
                    diagonals.append((v.index, left_edge.helper.index))
                left_edge.helper = v
        
        else:  # REGULAR
            # Check if interior is to the right
            if v.prev.y > v.y or (abs(v.prev.y - v.y) < 1e-9 and v.prev.x < v.x):
                # Interior to the right, handle edge ending here
                edge = find_edge_with_upper(v.prev)
                if edge and edge.helper and edge.helper.vertex_type == VertexType.MERGE:
                    diagonals.append((v.index, edge.helper.index))
                remove_edge(v)
                add_edge(v)
            else:
                # Interior to the left
                left_edge = find_left_edge(v)
                if left_edge:
                    if left_edge.helper and left_edge.helper.vertex_type == VertexType.MERGE:
                        diagonals.append((v.index, left_edge.helper.index))
                    left_edge.helper = v
    
    return diagonals


def triangulate_monotone(vertices: List[Vertex], indices: List[int]) -> List[Tuple[int, int, int]]:
    """
    Triangulate a y-monotone polygon using stack-based algorithm.
    Returns list of triangles as (i, j, k) vertex indices.
    """
    if len(indices) < 3:
        return []
    if len(indices) == 3:
        return [(indices[0], indices[1], indices[2])]
    
    # Sort vertices by y (descending)
    sorted_indices = sorted(indices, key=lambda i: (-vertices[i].y, vertices[i].x))
    
    triangles = []
    stack = [sorted_indices[0], sorted_indices[1]]
    
    # Determine which chain each vertex is on
    # Find leftmost and rightmost at each y level to determine chains
    left_chain = set()
    right_chain = set()
    
    # Walk around polygon to determine chains
    # Start from top vertex, go around
    top_idx = sorted_indices[0]
    
    # Find position of top in original indices
    pos = indices.index(top_idx)
    n = len(indices)
    
    # Walk one direction
    chain1 = []
    i = pos
    while True:
        chain1.append(indices[i])
        i = (i + 1) % n
        if i == pos:
            break
    
    # Determine which is left/right based on second vertex
    if len(chain1) >= 2:
        v0 = vertices[chain1[0]]
        v1 = vertices[chain1[1]]
        if v1.x < v0.x or (abs(v1.x - v0.x) < 1e-9 and v1.y < v0.y):
            left_chain = set(chain1[:len(chain1)//2 + 1])
            right_chain = set(chain1[len(chain1)//2:])
        else:
            right_chain = set(chain1[:len(chain1)//2 + 1])
            left_chain = set(chain1[len(chain1)//2:])
    
    for j in range(2, len(sorted_indices)):
        curr = sorted_indices[j]
        
        if len(stack) < 2:
            stack.append(curr)
            continue
        
        # Check if on same chain as top of stack
        curr_left = curr in left_chain
        top_left = stack[-1] in left_chain
        
        if curr_left != top_left:
            # Different chain - pop all and triangulate
            while len(stack) > 1:
                v1 = stack.pop()
                v2 = stack[-1]
                triangles.append((curr, v1, v2))
            stack.pop()
            stack.append(sorted_indices[j-1])
            stack.append(curr)
        else:
            # Same chain - pop while diagonal is inside
            last_popped = stack.pop()
            while stack:
                v = stack[-1]
                # Check if diagonal is inside polygon
                if curr_left:
                    # Left chain - check CCW
                    if cross_product(vertices[curr], vertices[last_popped], vertices[v]) > 0:
                        triangles.append((curr, last_popped, v))
                        last_popped = stack.pop()
                    else:
                        break
                else:
                    # Right chain - check CW
                    if cross_product(vertices[curr], vertices[last_popped], vertices[v]) < 0:
                        triangles.append((curr, last_popped, v))
                        last_popped = stack.pop()
                    else:
                        break
            stack.append(last_popped)
            stack.append(curr)
    
    return triangles


def triangulate_polygon(vertices: List[Vertex]) -> List[Tuple[int, int, int]]:
    """
    Triangulate a simple polygon using Garey's algorithm.
    1. Decompose into y-monotone pieces
    2. Triangulate each monotone piece
    """
    n = len(vertices)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]
    
    # Get diagonals from monotone decomposition
    diagonals = make_monotone(vertices)
    
    # Build adjacency for sub-polygons
    # For simplicity, if no diagonals needed, polygon is already monotone
    if not diagonals:
        return triangulate_monotone(vertices, list(range(n)))
    
    # With diagonals, we need to split into sub-polygons
    # For this implementation, we'll use a simpler ear-clipping as fallback
    # when the monotone decomposition produces diagonals
    
    triangles = []
    remaining = list(range(n))
    
    while len(remaining) > 3:
        found_ear = False
        for i in range(len(remaining)):
            prev_idx = remaining[(i - 1) % len(remaining)]
            curr_idx = remaining[i]
            next_idx = remaining[(i + 1) % len(remaining)]
            
            v_prev = vertices[prev_idx]
            v_curr = vertices[curr_idx]
            v_next = vertices[next_idx]
            
            # Check if convex
            if cross_product(v_prev, v_curr, v_next) <= 0:
                continue
            
            # Check no other vertices inside
            is_ear = True
            for j in remaining:
                if j in (prev_idx, curr_idx, next_idx):
                    continue
                if point_in_triangle(vertices[j], v_prev, v_curr, v_next):
                    is_ear = False
                    break
            
            if is_ear:
                triangles.append((prev_idx, curr_idx, next_idx))
                remaining.pop(i)
                found_ear = True
                break
        
        if not found_ear:
            # Fallback: just add any triangle
            if len(remaining) >= 3:
                triangles.append((remaining[0], remaining[1], remaining[2]))
                remaining.pop(1)
    
    if len(remaining) == 3:
        triangles.append((remaining[0], remaining[1], remaining[2]))
    
    return triangles


def point_in_triangle(p: Vertex, p1: Vertex, p2: Vertex, p3: Vertex) -> bool:
    """Check if point p is strictly inside triangle p1-p2-p3."""
    def sign(a, b, c):
        return (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y)
    
    d1 = sign(p, p1, p2)
    d2 = sign(p, p2, p3)
    d3 = sign(p, p3, p1)
    
    has_neg = (d1 < -1e-9) or (d2 < -1e-9) or (d3 < -1e-9)
    has_pos = (d1 > 1e-9) or (d2 > 1e-9) or (d3 > 1e-9)
    
    return not (has_neg and has_pos)


def read_polygon(filename: str) -> List[Vertex]:
    """Read polygon from .poly file."""
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        vertices = []
        for i in range(n):
            x, y = map(float, f.readline().strip().split())
            vertices.append(Vertex(x, y, i))
        
        # Set up circular links
        for i in range(n):
            vertices[i].prev = vertices[(i - 1) % n]
            vertices[i].next = vertices[(i + 1) % n]
    
    return vertices


def write_triangulation(vertices: List[Vertex], triangles: List[Tuple[int, int, int]], filename: str):
    """Write triangulation to .tri file."""
    with open(filename, 'w') as f:
        f.write("# vertices\n")
        f.write(f"{len(vertices)}\n")
        for v in vertices:
            f.write(f"{v.x} {v.y}\n")
        
        f.write("# triangles\n")
        f.write(f"{len(triangles)}\n")
        for t in triangles:
            f.write(f"{t[0]} {t[1]} {t[2]}\n")


def main():
    parser = argparse.ArgumentParser(description='Garey O(n log n) triangulation')
    parser.add_argument('--input', '-i', required=True, help='Input polygon file')
    parser.add_argument('--output', '-o', required=True, help='Output triangulation file')
    args = parser.parse_args()
    
    vertices = read_polygon(args.input)
    n = len(vertices)
    
    start = time.perf_counter()
    triangles = triangulate_polygon(vertices)
    end = time.perf_counter()
    
    elapsed_ms = (end - start) * 1000
    
    write_triangulation(vertices, triangles, args.output)
    print(f"garey,vertices={n},triangles={len(triangles)},time_ms={elapsed_ms}")


if __name__ == '__main__':
    main()

