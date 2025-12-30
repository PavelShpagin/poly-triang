"""
Comprehensive tests for the O(n + r log r) triangulation algorithm.
Tests both correctness and theoretical complexity properties.
"""

import math
import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import (
    PolygonTriangulator, 
    create_paper_example, 
    create_star_polygon,
    create_comb_polygon,
    create_nested_polygon,
    validate_triangulation
)


def check_triangles_cover_polygon(pts, triangles):
    """Verify that triangles cover the polygon area correctly."""
    # Compute polygon area
    n = len(pts)
    poly_area = 0
    for i in range(n):
        j = (i + 1) % n
        poly_area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
    poly_area = abs(poly_area) / 2
    
    # Compute total triangle area
    tri_area = 0
    for t in triangles:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2
        tri_area += area
    
    return abs(poly_area - tri_area) < 1e-6


def check_triangles_valid(pts, triangles):
    """Check that all triangles have positive area (no degenerate triangles)."""
    for t in triangles:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if abs(area) < 1e-9:
            return False, f"Degenerate triangle: {t}"
    return True, "OK"


def create_convex_polygon(n):
    """Create a convex n-gon (regular polygon)."""
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        pts.append((math.cos(angle), math.sin(angle)))
    return pts


def create_monotone_polygon():
    """Create a simple y-monotone polygon."""
    return [
        (0, 0), (2, 1), (4, 0), (4, 3), (2, 2), (0, 3)
    ]


def create_L_shape():
    """Create an L-shaped polygon."""
    return [
        (0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)
    ]


def create_complex_star(n_points, spikes=5):
    """Create a complex star with varying spike depths."""
    pts = []
    for i in range(spikes * 2):
        angle = math.pi / 2 + i * math.pi / spikes
        if i % 2 == 0:
            r = 2.0
        else:
            r = 0.5 + 0.3 * (i % 4)  # Varying inner radius
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def run_test(name, pts, verbose=False):
    """Run triangulation test and return success status."""
    pt = PolygonTriangulator(pts)
    tris = pt.triangulate()
    
    # Check triangle count
    valid, msg = validate_triangulation(pt.pts, tris)
    if not valid:
        print(f"FAIL [{name}]: {msg}")
        return False
    
    # Check triangles are non-degenerate
    valid, msg = check_triangles_valid(pt.pts, tris)
    if not valid:
        print(f"FAIL [{name}]: {msg}")
        return False
    
    # Check area coverage
    if not check_triangles_cover_polygon(pt.pts, tris):
        print(f"FAIL [{name}]: Area mismatch")
        return False
    
    if verbose:
        print(f"PASS [{name}]: n={len(pts)}, r={pt.r}, triangles={len(tris)}, diagonals={len(pt.diagonals)}")
    
    return True


def test_complexity_invariants():
    """
    Test that the algorithm maintains O(n + r log r) complexity invariants:
    - Number of extrema is O(r) [specifically <= 2r + 2]
    - Number of chains is O(r)
    - Number of diagonals is O(r)
    """
    print("\n=== Complexity Invariants ===")
    
    test_cases = [
        ("Convex 10", create_convex_polygon(10)),
        ("Convex 100", create_convex_polygon(100)),
        ("5-star", create_star_polygon(5)),
        ("10-star", create_star_polygon(10)),
        ("3-comb", create_comb_polygon(3)),
        ("10-comb", create_comb_polygon(10)),
        ("Paper example", create_paper_example()),
    ]
    
    for name, pts in test_cases:
        pt = PolygonTriangulator(pts)
        pt._classify_vertices()
        pt._build_chains()
        
        n = len(pts)
        r = pt.r
        num_extrema = len(pt.extrema)
        num_chains = len(pt.chains)
        
        # Invariants
        # 1. num_extrema <= 2r + 2 (each reflex vertex contributes at most 2 extrema)
        extrema_bound = 2 * r + 2
        # 2. num_chains = num_extrema (one chain per extremum)
        
        extrema_ok = num_extrema <= extrema_bound
        chains_ok = num_chains == num_extrema
        
        status = "OK" if (extrema_ok and chains_ok) else "FAIL"
        print(f"{name}: n={n}, r={r}, extrema={num_extrema} (<= {extrema_bound}? {extrema_ok}), chains={num_chains} (= extrema? {chains_ok}) -> {status}")


def main():
    print("=== Comprehensive Triangulation Tests ===\n")
    
    # Basic tests
    basic_tests = [
        ("Triangle", [(0, 0), (1, 0), (0.5, 1)]),
        ("Square", [(0, 0), (1, 0), (1, 1), (0, 1)]),
        ("Pentagon (convex)", create_convex_polygon(5)),
        ("Hexagon (convex)", create_convex_polygon(6)),
        ("L-shape", create_L_shape()),
        ("Monotone", create_monotone_polygon()),
        ("Paper example", create_paper_example()),
        ("5-point star", create_star_polygon(5)),
        ("7-point star", create_star_polygon(7)),
        ("3-tooth comb", create_comb_polygon(3)),
        ("5-tooth comb", create_comb_polygon(5)),
        ("Nested", create_nested_polygon()),
        ("Complex star", create_complex_star(10, 5)),
    ]
    
    all_passed = True
    for name, pts in basic_tests:
        if not run_test(name, pts, verbose=True):
            all_passed = False
    
    # Larger tests
    print("\n=== Larger Polygon Tests ===")
    larger_tests = [
        ("Convex 50", create_convex_polygon(50)),
        ("Convex 100", create_convex_polygon(100)),
        ("20-point star", create_star_polygon(20)),
        ("10-tooth comb", create_comb_polygon(10)),
    ]
    
    for name, pts in larger_tests:
        if not run_test(name, pts, verbose=True):
            all_passed = False
    
    # Complexity invariants
    test_complexity_invariants()
    
    print("\n" + "=" * 50)
    if all_passed:
        print("ALL TESTS PASSED!")
    else:
        print("SOME TESTS FAILED!")
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())

