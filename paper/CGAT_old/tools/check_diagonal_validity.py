#!/usr/bin/env python3
"""
Check validity of diagonals produced by chain_only on a generated polygon.

We look for:
 - diagonal intersecting polygon boundary (excluding incident edges)
 - diagonal intersecting another diagonal (excluding shared endpoints)
"""

import math
import random
import re
import subprocess
import os
import argparse
from typing import List, Tuple


ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def simple_random_poly(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    rng = random.Random(seed + n * 1000)
    angles = sorted([rng.random() * 2 * math.pi for _ in range(n)])
    pts = []
    for angle in angles:
        r = 50 + rng.random() * 50
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def write_poly(pts: List[Tuple[float, float]], path: str) -> None:
    with open(path, "w") as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")


EPS = 1e-12


def orient(ax, ay, bx, by, cx, cy) -> float:
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def on_segment(ax, ay, bx, by, cx, cy) -> bool:
    # c collinear with a-b assumed
    return (
        min(ax, bx) - 1e-10 <= cx <= max(ax, bx) + 1e-10
        and min(ay, by) - 1e-10 <= cy <= max(ay, by) + 1e-10
    )


def segments_properly_intersect(a, b, c, d) -> bool:
    """True if segments ab and cd intersect at a point that is not a shared endpoint."""
    ax, ay = a
    bx, by = b
    cx, cy = c
    dx, dy = d

    o1 = orient(ax, ay, bx, by, cx, cy)
    o2 = orient(ax, ay, bx, by, dx, dy)
    o3 = orient(cx, cy, dx, dy, ax, ay)
    o4 = orient(cx, cy, dx, dy, bx, by)

    def sgn(x):
        if x > EPS:
            return 1
        if x < -EPS:
            return -1
        return 0

    s1, s2, s3, s4 = sgn(o1), sgn(o2), sgn(o3), sgn(o4)

    # General proper intersection
    if s1 * s2 < 0 and s3 * s4 < 0:
        return True

    # Collinear / touching cases: treat as intersection if overlapping beyond endpoints
    if s1 == 0 and on_segment(ax, ay, bx, by, cx, cy):
        return True
    if s2 == 0 and on_segment(ax, ay, bx, by, dx, dy):
        return True
    if s3 == 0 and on_segment(cx, cy, dx, dy, ax, ay):
        return True
    if s4 == 0 and on_segment(cx, cy, dx, dy, bx, by):
        return True

    return False


def parse_diagonals_from_diag_debug(out: str) -> List[Tuple[int, int]]:
    diags = []
    in_diag = False
    in_raw = False
    for line in out.splitlines():
        if line.startswith("# diagonals:"):
            in_raw = True
            in_diag = False
            continue
        if line.startswith("# diag_records:"):
            in_diag = True
            in_raw = False
            continue
        if line.startswith("# pending_records:"):
            in_diag = False
            in_raw = False
        if not in_diag:
            if not in_raw:
                continue
        parts = line.strip().split(",")
        if len(parts) < 2:
            continue
        try:
            a = int(parts[0])
            b = int(parts[1])
        except Exception:
            continue
        diags.append((a, b))
    # unique
    diags = list({(min(a, b), max(a, b)) for a, b in diags})
    diags.sort()
    return diags


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--poly", type=str, default="")
    args = ap.parse_args()

    if args.poly:
        poly_path = args.poly
        with open(poly_path) as f:
            n = int(f.readline().strip())
            pts = [tuple(map(float, f.readline().split())) for _ in range(n)]
    else:
        n, seed = args.n, args.seed
        pts = simple_random_poly(n, seed)
        poly_path = "/tmp/diag_check.poly"
        write_poly(pts, poly_path)

    # Prefer the reproducible build output (paper/CGAT/build.sh builds this).
    diag_dbg = os.path.join(ROOT, "bin", "diag_debug_cli")
    if not os.path.exists(diag_dbg):
        # Fallback for local dev builds (not used in the artifact pipeline).
        diag_dbg = os.path.join(ROOT, "code", "ours", "diag_debug_cli")
    if not os.path.exists(diag_dbg):
        print(f"diag_debug_cli not found at '{os.path.join(ROOT, 'bin', 'diag_debug_cli')}'.")
        print("Run ./build.sh (or ./reproduce.sh) from paper/CGAT first.")
        return 2

    res = subprocess.run([diag_dbg, "--input", poly_path], capture_output=True, text=True)
    if res.returncode != 0:
        print("diag_debug_cli failed:")
        print(res.stdout[:500])
        print(res.stderr[:500])
        return 2
    diags = parse_diagonals_from_diag_debug(res.stdout)
    if not diags:
        print("Failed to parse diagonals from diag_debug_cli output.")
        print(res.stdout[:500])
        return 2

    # Boundary edges
    boundary = set()
    for i in range(n):
        j = (i + 1) % n
        boundary.add((min(i, j), max(i, j)))

    # Check diagonals against boundary intersections
    bad_boundary = []
    for (u, v) in diags:
        a = pts[u]
        b = pts[v]
        # Skip if diagonal is actually a boundary edge
        if (min(u, v), max(u, v)) in boundary:
            bad_boundary.append(((u, v), "is boundary edge"))
            continue

        for i in range(n):
            j = (i + 1) % n
            # Ignore edges incident to u or v
            if i in (u, v) or j in (u, v):
                continue
            if segments_properly_intersect(a, b, pts[i], pts[j]):
                bad_boundary.append(((u, v), f"intersects boundary edge ({i},{j})"))
                break

    # Check diagonal-diagonal crossings
    bad_cross = []
    for i in range(len(diags)):
        u1, v1 = diags[i]
        a1, b1 = pts[u1], pts[v1]
        for j in range(i + 1, len(diags)):
            u2, v2 = diags[j]
            # Shared endpoint allowed
            if len({u1, v1, u2, v2}) < 4:
                continue
            a2, b2 = pts[u2], pts[v2]
            if segments_properly_intersect(a1, b1, a2, b2):
                bad_cross.append(((u1, v1), (u2, v2)))

    if args.poly:
        print(f"poly={poly_path}, n={n}, diagonals={len(diags)}")
    else:
        print(f"n={n}, seed={args.seed}, diagonals={len(diags)}")
    print(f"boundary-issues: {len(bad_boundary)}")
    for d, why in bad_boundary[:20]:
        print(f"  diag {d}: {why}")
    if len(bad_boundary) > 20:
        print("  ...")

    print(f"diagonal-crossings: {len(bad_cross)}")
    for d1, d2 in bad_cross[:20]:
        print(f"  {d1} x {d2}")
    if len(bad_cross) > 20:
        print("  ...")

    # Treat any reported issue as a failure (artifact should fail-fast on correctness bugs).
    if bad_boundary or bad_cross:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

