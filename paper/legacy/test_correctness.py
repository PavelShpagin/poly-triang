#!/usr/bin/env python3
"""
Correctness verification for the O(n + r log r) triangulation algorithm.

Tests:
1. Triangle count: n - 2 triangles for n-vertex polygon
2. Area preservation: sum of triangle areas == polygon area
3. No degenerate triangles: all triangles have positive area
4. Valid indices: all triangle vertices are valid polygon indices
"""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple

ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = ROOT / "build" / "bin"


def polygon_area(pts: List[Tuple[float, float]]) -> float:
    """Compute signed area of polygon."""
    n = len(pts)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
    return abs(area) / 2


def triangle_area(pts: List[Tuple[float, float]], tri: Tuple[int, int, int]) -> float:
    """Compute signed area of triangle."""
    a, b, c = pts[tri[0]], pts[tri[1]], pts[tri[2]]
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    with open(path, "w") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            f.write(f"{x:.10f} {y:.10f}\n")


def read_triangles(path: Path) -> Tuple[List[Tuple[float, float]], List[Tuple[int, int, int]]]:
    """Read triangulation output."""
    pts = []
    tris = []
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    
    i = 0
    n = int(lines[i])
    i += 1
    for _ in range(n):
        x, y = map(float, lines[i].split())
        pts.append((x, y))
        i += 1
    
    if i < len(lines):
        nt = int(lines[i])
        i += 1
        for _ in range(nt):
            v0, v1, v2 = map(int, lines[i].split())
            tris.append((v0, v1, v2))
            i += 1
    
    return pts, tris


def run_triangulation(exe: Path, pts: List[Tuple[float, float]]) -> Tuple[bool, List[Tuple[int, int, int]], str]:
    """Run triangulation and return (success, triangles, error_msg)."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        poly_path = td / "poly.poly"
        out_path = td / "out.tri"
        
        write_poly(pts, poly_path)
        
        try:
            proc = subprocess.run(
                [str(exe), "--input", str(poly_path), "--output", str(out_path)],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if proc.returncode != 0:
                return False, [], proc.stderr[:200]
            
            _, tris = read_triangles(out_path)
            return True, tris, ""
        except Exception as e:
            return False, [], str(e)


def verify_triangulation(pts: List[Tuple[float, float]], tris: List[Tuple[int, int, int]]) -> Tuple[bool, str]:
    """Verify that triangulation is correct."""
    n = len(pts)
    
    # Check triangle count
    expected = n - 2
    if len(tris) != expected:
        return False, f"Wrong count: {len(tris)} != {expected}"
    
    # Check valid indices
    for tri in tris:
        for v in tri:
            if v < 0 or v >= n:
                return False, f"Invalid vertex index: {v}"
    
    # Check no degenerate triangles
    for tri in tris:
        area = triangle_area(pts, tri)
        if area < 1e-12:
            return False, f"Degenerate triangle: {tri}"
    
    # Check area preservation
    poly_a = polygon_area(pts)
    tri_a = sum(triangle_area(pts, tri) for tri in tris)
    if abs(poly_a - tri_a) > 1e-6 * max(1, poly_a):
        return False, f"Area mismatch: {poly_a:.6f} vs {tri_a:.6f}"
    
    return True, "OK"


# Polygon generators
def convex_polygon(n: int) -> List[Tuple[float, float]]:
    return [(math.cos(2 * math.pi * i / n), math.sin(2 * math.pi * i / n)) for i in range(n)]


def star_polygon(points: int) -> List[Tuple[float, float]]:
    pts = []
    for i in range(points * 2):
        angle = math.pi / 2 + i * math.pi / points
        r = 2.0 if i % 2 == 0 else 0.8
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def l_shape() -> List[Tuple[float, float]]:
    return [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]


def arrow_shape() -> List[Tuple[float, float]]:
    return [(0, 1), (2, 1), (2, 0), (4, 1.5), (2, 3), (2, 2), (0, 2)]


def comb_polygon(teeth: int) -> List[Tuple[float, float]]:
    pts = [(0, 0), (teeth * 2, 0), (teeth * 2, 1)]
    for i in range(teeth - 1, -1, -1):
        x = i * 2 + 1
        pts.extend([(x + 0.5, 1), (x, 2), (x - 0.5, 1)])
    pts.append((0, 1))
    return pts


def paper_example() -> List[Tuple[float, float]]:
    """The example from the paper."""
    return [
        (0.0, 2.5), (1.2, 5.5), (2.5, 3.8), (4.0, 6.5),
        (5.5, 4.8), (7.0, 7.0), (8.0, 5.5), (6.5, 3.5),
        (8.0, 1.5), (5.0, 2.5), (3.0, 0.0), (1.5, 1.5),
    ]


def spiral_polygon(n: int) -> List[Tuple[float, float]]:
    """Create a simple (non-self-intersecting) spiral-like polygon.
    
    This creates a star-burst pattern with n/2 outer points and n/2 inner points,
    which is always simple and has many reflex vertices.
    """
    if n < 6:
        n = 6
    half = n // 2
    pts = []
    for i in range(half):
        angle = 2 * math.pi * i / half
        # Outer point
        pts.append((100 * math.cos(angle), 100 * math.sin(angle)))
        # Inner point (offset by half step)
        inner_angle = angle + math.pi / half
        pts.append((30 * math.cos(inner_angle), 30 * math.sin(inner_angle)))
    return pts


def random_polygon(n: int, seed: int = 42) -> List[Tuple[float, float]]:
    # Use the repo's dataset generator model (same as scripts/generate_polygons.py),
    # which is intended to produce simple polygons in general position.
    import random
    rng = random.Random(seed + n)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    pts = []
    for a in angles:
        r = 100.0 * (0.4 + 0.6 * rng.random())
        pts.append((r * math.cos(a), r * math.sin(a)))
    # Deterministic tiny rotation to avoid equal-y ties after rounding in C++.
    rot = 0.123456789
    ca, sa = math.cos(rot), math.sin(rot)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in pts]


def main():
    print("=" * 60)
    print("Triangulation Correctness Verification")
    print("=" * 60)
    
    exe = BIN_DIR / "reflex_cli"
    if not exe.exists():
        print(f"ERROR: {exe} not found. Build first.")
        return 1
    
    test_cases = [
        ("Triangle", [(0, 0), (1, 0), (0.5, 1)]),
        ("Square", [(0, 0), (1, 0), (1, 1), (0, 1)]),
        ("Pentagon (convex)", convex_polygon(5)),
        ("Hexagon (convex)", convex_polygon(6)),
        ("Decagon (convex)", convex_polygon(10)),
        ("L-shape", l_shape()),
        ("Arrow", arrow_shape()),
        ("Paper example", paper_example()),
        ("5-star", star_polygon(5)),
        ("6-star", star_polygon(6)),
        ("7-star", star_polygon(7)),
        ("10-star", star_polygon(10)),
        ("Comb-3", comb_polygon(3)),
        ("Comb-5", comb_polygon(5)),
        ("Comb-10", comb_polygon(10)),
        ("Convex 50", convex_polygon(50)),
        ("Convex 100", convex_polygon(100)),
        ("Random 50", random_polygon(50)),
        ("Random 100", random_polygon(100)),
        ("Random 500", random_polygon(500)),
        ("Spiral 50", spiral_polygon(50)),
        ("Spiral 100", spiral_polygon(100)),
    ]
    
    passed = 0
    failed = 0
    
    for name, pts in test_cases:
        success, tris, error = run_triangulation(exe, pts)
        
        if not success:
            print(f"FAIL [{name}]: Execution error - {error}")
            failed += 1
            continue
        
        valid, msg = verify_triangulation(pts, tris)
        
        if valid:
            print(f"PASS [{name}]: n={len(pts)}, triangles={len(tris)}")
            passed += 1
        else:
            print(f"FAIL [{name}]: {msg}")
            failed += 1
    
    print()
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

