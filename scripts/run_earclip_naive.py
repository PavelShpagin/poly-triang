#!/usr/bin/env python3
"""
Naive O(n^2) Ear Clipping implementation in Python.
No optimizations - pure textbook algorithm to show true O(n^2) behavior.
"""

import sys
import time
import argparse
from typing import List, Tuple


class Vertex:
    def __init__(self, x: float, y: float, index: int):
        self.x = x
        self.y = y
        self.index = index
        self.prev = None
        self.next = None


def cross_product(o: Vertex, a: Vertex, b: Vertex) -> float:
    """Cross product of vectors OA and OB."""
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def point_in_triangle(p: Vertex, a: Vertex, b: Vertex, c: Vertex) -> bool:
    """Check if point p is strictly inside triangle abc."""
    def sign(p1, p2, p3):
        return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
    
    d1 = sign(p, a, b)
    d2 = sign(p, b, c)
    d3 = sign(p, c, a)
    
    has_neg = (d1 < -1e-10) or (d2 < -1e-10) or (d3 < -1e-10)
    has_pos = (d1 > 1e-10) or (d2 > 1e-10) or (d3 > 1e-10)
    
    return not (has_neg and has_pos)


def is_convex(v: Vertex) -> bool:
    """Check if vertex v is convex (interior angle < 180 degrees)."""
    return cross_product(v.prev, v, v.next) > 0


def is_ear(v: Vertex, vertices: List[Vertex]) -> bool:
    """
    Check if vertex v is an ear.
    An ear is a convex vertex where no other polygon vertices are inside the triangle.
    This is the O(n) check that makes the algorithm O(n^2) overall.
    """
    if not is_convex(v):
        return False
    
    a, b, c = v.prev, v, v.next
    
    # Check all other vertices - O(n) operation
    for other in vertices:
        if other is a or other is b or other is c:
            continue
        if other.prev is None:  # Already removed
            continue
        if point_in_triangle(other, a, b, c):
            return False
    
    return True


def ear_clip_naive(vertices: List[Vertex]) -> List[Tuple[int, int, int]]:
    """
    Naive ear clipping triangulation - O(n^2) algorithm.
    
    For each of the n-2 triangles:
      - Scan all remaining vertices to find an ear: O(n)
      - For each candidate, check all vertices for point-in-triangle: O(n)
    Total: O(n^2) worst case
    """
    triangles = []
    n = len(vertices)
    remaining = n
    
    # Set up doubly linked list
    for i in range(n):
        vertices[i].prev = vertices[(i - 1) % n]
        vertices[i].next = vertices[(i + 1) % n]
    
    # Need to produce n-2 triangles
    current = vertices[0]
    iterations = 0
    max_iterations = n * n  # Safety limit
    
    while remaining > 3 and iterations < max_iterations:
        iterations += 1
        found_ear = False
        
        # Scan for an ear - start from current vertex
        start = current
        v = current
        
        while True:
            if is_ear(v, vertices):
                # Found an ear - add triangle
                triangles.append((v.prev.index, v.index, v.next.index))
                
                # Remove vertex from linked list
                v.prev.next = v.next
                v.next.prev = v.prev
                v.prev = None  # Mark as removed
                
                remaining -= 1
                current = v.next
                found_ear = True
                break
            
            v = v.next
            if v is start:
                break
        
        if not found_ear:
            # No ear found - should not happen for valid simple polygon
            # Try next vertex anyway
            current = current.next
    
    # Add final triangle
    if remaining == 3:
        v = current
        triangles.append((v.prev.index, v.index, v.next.index))
    
    return triangles


def read_polygon(filename: str) -> List[Vertex]:
    """Read polygon from .poly file."""
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        vertices = []
        for i in range(n):
            x, y = map(float, f.readline().strip().split())
            vertices.append(Vertex(x, y, i))
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
    parser = argparse.ArgumentParser(description='Naive O(n^2) Ear Clipping')
    parser.add_argument('--input', '-i', required=True, help='Input polygon file')
    parser.add_argument('--output', '-o', required=True, help='Output triangulation file')
    args = parser.parse_args()
    
    vertices = read_polygon(args.input)
    n = len(vertices)
    
    start = time.perf_counter()
    triangles = ear_clip_naive(vertices)
    end = time.perf_counter()
    
    elapsed_ms = (end - start) * 1000
    
    write_triangulation(vertices, triangles, args.output)
    print(f"earclip_naive,vertices={n},triangles={len(triangles)},time_ms={elapsed_ms}")


if __name__ == '__main__':
    main()

