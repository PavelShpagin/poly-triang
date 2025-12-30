#!/usr/bin/env python3
"""
Generate deterministic polygon datasets for benchmarking.
The format is:
N
x0 y0
x1 y1
...
"""

import argparse
import math
import random
from pathlib import Path


def write_polygon(points, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            # Use high precision to avoid accidental equal y after rounding.
            f.write(f"{x:.17g} {y:.17g}\n")


def rotate_points(points, angle_rad):
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in points]


def convex_polygon(n, radius=100.0):
    return [
        (
            radius * math.cos(2 * math.pi * i / n),
            radius * math.sin(2 * math.pi * i / n),
        )
        for i in range(n)
    ]


def random_polygon(n, radius=100.0, seed=42):
    rng = random.Random(seed + n)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    points = []
    for angle in angles:
        r = radius * (0.4 + 0.6 * rng.random())
        points.append((r * math.cos(angle), r * math.sin(angle)))
    return points


def star_polygon(n_pairs, outer=100.0, inner=30.0):
    points = []
    for i in range(2 * n_pairs):
        angle = math.pi * i / n_pairs
        r = outer if i % 2 == 0 else inner
        points.append((r * math.cos(angle), r * math.sin(angle)))
    return points


def spiral_polygon(n, turns=3, start=20.0, end=100.0):
    points = []
    for i in range(n):
        t = i / (n - 1)
        angle = 2 * math.pi * turns * t
        radius = start + (end - start) * t
        points.append((radius * math.cos(angle), radius * math.sin(angle)))
    return points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="polygons/generated", type=Path)
    parser.add_argument(
        "--sizes",
        nargs="+",
        type=int,
        default=[10, 50, 100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000],
    )
    args = parser.parse_args()

    # Deterministic rotation to avoid general-position degeneracies (equal y) in saved files.
    # The paper algorithm assumes general position; this keeps datasets consistent with it.
    rot = 0.123456789  # radians

    for n in args.sizes:
        write_polygon(rotate_points(convex_polygon(n), rot), args.output / f"convex_{n}.poly")
        write_polygon(rotate_points(random_polygon(n), rot), args.output / f"random_{n}.poly")
        write_polygon(rotate_points(star_polygon(max(3, n // 2)), rot), args.output / f"star_{n}.poly")
        write_polygon(rotate_points(spiral_polygon(n), rot), args.output / f"spiral_{n}.poly")

    print(f"Generated polygons in {args.output}")


if __name__ == "__main__":
    main()

