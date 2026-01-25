#!/usr/bin/env python3
"""
Validate a triangulation in CGAT .tri format against the polygon boundary.

Checks:
 - triangles == n-2
 - all indices in range, no degenerate triangles
 - no edge crossings (triangle edges vs triangle edges, excluding shared endpoints)
 - no triangle edge (except boundary edges) intersects polygon boundary (excluding shared endpoints)

This is intended for *reviewer confidence* on the small visualization instances (n≈200),
not as part of the high-performance benchmark loop.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple, Set

EPS = 1e-12


def read_tri(path: Path) -> Tuple[List[Tuple[float, float]], List[Tuple[int, int, int]]]:
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines or lines[0] != "# vertices":
        raise RuntimeError(f"Unexpected .tri format (missing '# vertices'): {path}")
    n = int(lines[1])
    pts = [tuple(map(float, lines[2 + i].split())) for i in range(n)]
    k = 2 + n
    if k >= len(lines) or lines[k] != "# triangles":
        raise RuntimeError(f"Unexpected .tri format (missing '# triangles'): {path}")
    m = int(lines[k + 1])
    tris = [tuple(map(int, lines[k + 2 + i].split())) for i in range(m)]
    return pts, tris


def orient(a: Tuple[float, float], b: Tuple[float, float], c: Tuple[float, float]) -> float:
    ax, ay = a
    bx, by = b
    cx, cy = c
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def on_segment(a: Tuple[float, float], b: Tuple[float, float], c: Tuple[float, float]) -> bool:
    ax, ay = a
    bx, by = b
    cx, cy = c
    return (
        min(ax, bx) - 1e-10 <= cx <= max(ax, bx) + 1e-10
        and min(ay, by) - 1e-10 <= cy <= max(ay, by) + 1e-10
    )


def segments_intersect(a: Tuple[float, float], b: Tuple[float, float],
                       c: Tuple[float, float], d: Tuple[float, float]) -> bool:
    """Counts touching/collinear overlap as intersection."""
    o1 = orient(a, b, c)
    o2 = orient(a, b, d)
    o3 = orient(c, d, a)
    o4 = orient(c, d, b)

    def sgn(x: float) -> int:
        if x > EPS:
            return 1
        if x < -EPS:
            return -1
        return 0

    s1, s2, s3, s4 = sgn(o1), sgn(o2), sgn(o3), sgn(o4)
    if s1 * s2 < 0 and s3 * s4 < 0:
        return True
    if s1 == 0 and on_segment(a, b, c):
        return True
    if s2 == 0 and on_segment(a, b, d):
        return True
    if s3 == 0 and on_segment(c, d, a):
        return True
    if s4 == 0 and on_segment(c, d, b):
        return True
    return False


def undirected(u: int, v: int) -> Tuple[int, int]:
    return (u, v) if u < v else (v, u)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tri", type=Path, required=True)
    args = ap.parse_args()

    pts, tris = read_tri(args.tri)
    n = len(pts)

    ok = True
    if len(tris) != n - 2:
        print(f"FAIL: triangle count {len(tris)} != n-2 ({n-2})")
        ok = False

    boundary: Set[Tuple[int, int]] = set()
    for i in range(n):
        boundary.add(undirected(i, (i + 1) % n))

    edges: Set[Tuple[int, int]] = set()
    for t_idx, (a, b, c) in enumerate(tris):
        if not (0 <= a < n and 0 <= b < n and 0 <= c < n):
            print(f"FAIL: triangle {t_idx} has out-of-range index: {(a,b,c)}")
            ok = False
            continue
        if len({a, b, c}) != 3:
            print(f"FAIL: triangle {t_idx} has duplicate vertices: {(a,b,c)}")
            ok = False
            continue
        edges.add(undirected(a, b))
        edges.add(undirected(b, c))
        edges.add(undirected(c, a))

    edges_list = sorted(edges)

    # edge-edge crossings
    for i, (u1, v1) in enumerate(edges_list):
        a = pts[u1]
        b = pts[v1]
        for j in range(i + 1, len(edges_list)):
            u2, v2 = edges_list[j]
            if len({u1, v1, u2, v2}) < 4:
                continue
            c = pts[u2]
            d = pts[v2]
            if segments_intersect(a, b, c, d):
                print(f"FAIL: edge crossing {(u1,v1)} x {(u2,v2)}")
                return 1

    # non-boundary edges vs boundary
    boundary_list = [(i, (i + 1) % n) for i in range(n)]
    for (u, v) in edges_list:
        if (u, v) in boundary:
            continue
        a = pts[u]
        b = pts[v]
        for i, j in boundary_list:
            if i in (u, v) or j in (u, v):
                continue
            if segments_intersect(a, b, pts[i], pts[j]):
                print(f"FAIL: interior edge {(u,v)} intersects boundary edge {(i,j)}")
                return 1

    if ok:
        print(f"OK: n={n}, triangles={len(tris)}, edges={len(edges_list)}")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

