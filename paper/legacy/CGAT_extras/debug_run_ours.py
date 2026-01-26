#!/usr/bin/env python3
"""
Debug helper: generate one polygon instance and run CGAT `reflex_cli`.

This is intentionally small and dependency-free.
"""

from __future__ import annotations

import argparse
import math
import random
import subprocess
from pathlib import Path
from typing import List, Tuple

ROT_ANGLE = 0.123456789


def rotate_points(points: List[Tuple[float, float]], angle_rad: float) -> List[Tuple[float, float]]:
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in points]


def convex_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    r = 100.0
    pts = [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]
    pts = rotate_points(pts, ROT_ANGLE)
    rng = random.Random(seed + 1337 * n)
    sx = 1.0 + 0.05 * (2.0 * rng.random() - 1.0)
    sy = 1.0 + 0.05 * (2.0 * rng.random() - 1.0)
    shx = 0.05 * (2.0 * rng.random() - 1.0)
    shy = 0.02 * (2.0 * rng.random() - 1.0)
    ang = 0.15 * (2.0 * rng.random() - 1.0)
    pts = rotate_points(pts, ang)
    return [(sx * x + shx * y, shy * x + sy * y) for (x, y) in pts]


def dent_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    pts = convex_polygon(n, seed)
    if n < 5:
        return pts
    imax = max(range(n), key=lambda i: (pts[i][1], -pts[i][0], -i))
    cx = sum(x for x, _ in pts) / n
    cy = sum(y for _, y in pts) / n
    rng = random.Random(seed + 4242 * n)
    depth = 0.55 + 0.25 * rng.random()

    def apply_depth(d: float) -> Tuple[float, float]:
        x, y = pts[imax]
        return (x + d * (cx - x), y + d * (cy - y))

    p = (imax - 1 + n) % n
    nx = (imax + 1) % n
    new_pt = apply_depth(depth)
    safe = 0
    while safe < 20 and not (new_pt[1] < pts[p][1] and new_pt[1] < pts[nx][1]):
        depth = min(0.95, depth + 0.05)
        new_pt = apply_depth(depth)
        safe += 1
    pts[imax] = new_pt
    return pts


def random_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    rng = random.Random(seed + n)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    points: List[Tuple[float, float]] = []
    for angle in angles:
        rr = 100.0 * (0.4 + 0.6 * rng.random())
        points.append((rr * math.cos(angle), rr * math.sin(angle)))
    return rotate_points(points, ROT_ANGLE)


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            f.write(f"{x:.17g} {y:.17g}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--type", choices=["convex", "dent", "random"], default="random")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bin", type=Path, default=Path(__file__).resolve().parent.parent / "bin" / "reflex_cli")
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/cgat_debug"))
    args = ap.parse_args()

    gens = {"convex": convex_polygon, "dent": dent_polygon, "random": random_polygon}
    pts = gens[args.type](args.n, args.seed)
    poly = args.outdir / f"{args.type}_{args.n}_{args.seed}.poly"
    tri = args.outdir / f"{args.type}_{args.n}_{args.seed}.tri"
    write_poly(pts, poly)

    proc = subprocess.run(
        [str(args.bin), "--input", str(poly), "--output", str(tri), "--algo", "chain_only"],
        capture_output=True,
        text=True,
    )
    print("rc:", proc.returncode)
    if proc.stdout:
        print("stdout:", proc.stdout.strip()[:500])
    if proc.stderr:
        print("stderr:", proc.stderr.strip()[:2000])
    print("poly:", poly)
    print("tri:", tri if tri.exists() else "(not written)")


if __name__ == "__main__":
    main()

