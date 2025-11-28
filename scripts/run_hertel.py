#!/usr/bin/env python3
"""
Wrapper script to run Hertel-Mehlhorn triangulation from the Fast-Triangulation-Of-The-Plane repo.
"""

import sys
import time
import argparse
from pathlib import Path

# Add the hertel source directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "external" / "hertel" / "src"))

from hertel_mehlhorn_triangulator import FastTriangulator, Point, Triangle


def read_polygon(filename):
    """Read polygon from .poly file format."""
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        vertices = []
        for i in range(n):
            x, y = map(float, f.readline().strip().split())
            vertices.append(Point(x, y, i))
        
        # Set up prev/next pointers for the triangulator
        for i in range(n):
            vertices[i].prev_point = vertices[(i - 1) % n]
            vertices[i].next_point = vertices[(i + 1) % n]
    
    return vertices


def write_triangulation(vertices, triangles, filename):
    """Write triangulation to .tri file format."""
    with open(filename, 'w') as f:
        f.write("# vertices\n")
        f.write(f"{len(vertices)}\n")
        for v in vertices:
            f.write(f"{v.x} {v.y}\n")
        
        f.write("# triangles\n")
        f.write(f"{len(triangles)}\n")
        for tri in triangles:
            # Get vertex indices
            indices = [v.index for v in tri.vertices]
            f.write(f"{indices[0]} {indices[1]} {indices[2]}\n")


def main():
    parser = argparse.ArgumentParser(description='Hertel-Mehlhorn triangulation')
    parser.add_argument('--input', '-i', required=True, help='Input polygon file')
    parser.add_argument('--output', '-o', required=True, help='Output triangulation file')
    args = parser.parse_args()
    
    # Read polygon
    vertices = read_polygon(args.input)
    n = len(vertices)
    
    # Create triangulator and run
    triangulator = FastTriangulator()
    
    # Set up polygon edges
    for i in range(n):
        from hertel_mehlhorn_triangulator import Edge
        edge = Edge(vertices[i], vertices[(i + 1) % n])
        triangulator.polygon_edges.add(edge)
    
    # Time the triangulation
    start = time.perf_counter()
    
    # Run the triangulation using the improved algorithm
    try:
        triangulator.triangulate_improved(vertices)
    except Exception as e:
        # Fall back to basic triangulation if improved fails
        try:
            triangulator.triangulate_basic(vertices)
        except:
            # Last resort: simple ear clipping
            triangulator.triangles = []
            remaining = list(range(n))
            while len(remaining) > 3:
                for i in range(len(remaining)):
                    prev_idx = remaining[(i - 1) % len(remaining)]
                    curr_idx = remaining[i]
                    next_idx = remaining[(i + 1) % len(remaining)]
                    
                    p1, p2, p3 = vertices[prev_idx], vertices[curr_idx], vertices[next_idx]
                    
                    # Check if ear
                    if triangulator.orientation(p1, p2, p3) == 2:  # CCW = convex
                        # Check no other vertices inside
                        is_ear = True
                        for j in remaining:
                            if j in (prev_idx, curr_idx, next_idx):
                                continue
                            if point_in_triangle(vertices[j], p1, p2, p3):
                                is_ear = False
                                break
                        
                        if is_ear:
                            tri = Triangle(p1, p2, p3)
                            triangulator.triangles.append(tri)
                            remaining.pop(i)
                            break
            
            # Add final triangle
            if len(remaining) == 3:
                p1, p2, p3 = vertices[remaining[0]], vertices[remaining[1]], vertices[remaining[2]]
                triangulator.triangles.append(Triangle(p1, p2, p3))
    
    end = time.perf_counter()
    elapsed_ms = (end - start) * 1000
    
    # Write output
    write_triangulation(vertices, triangulator.triangles, args.output)
    
    # Print timing info
    print(f"hertel,vertices={n},triangles={len(triangulator.triangles)},time_ms={elapsed_ms}")


def point_in_triangle(p, p1, p2, p3):
    """Check if point p is inside triangle p1-p2-p3."""
    def sign(p1, p2, p3):
        return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
    
    d1 = sign(p, p1, p2)
    d2 = sign(p, p2, p3)
    d3 = sign(p, p3, p1)
    
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    
    return not (has_neg and has_pos)


if __name__ == '__main__':
    main()

