#!/usr/bin/env python3
"""
Generate simple SVG visualizations of triangulations produced by the CGAT binaries.

No third-party dependencies: writes plain SVG (polygon boundary + triangle edges).

Usage (typical from paper/CGAT):
  python3 tools/visualize_triangulation.py --bin ./bin/reflex_cli --outdir generated/figures

Outputs:
  - triangulation_convex.svg
  - triangulation_dent.svg
  - triangulation_random.svg
"""

from __future__ import annotations

import argparse
import math
import random
import subprocess
from dataclasses import dataclass
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


@dataclass(frozen=True)
class Triangulation:
    verts: List[Tuple[float, float]]
    tris: List[Tuple[int, int, int]]


def read_tri(path: Path) -> Triangulation:
    # Format written by CGAT CLIs:
    #   # vertices
    #   n
    #   x y  (n lines)
    #   # triangles
    #   m
    #   i j k (m lines)
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines or lines[0] != "# vertices":
        raise RuntimeError(f"Unexpected .tri format (missing '# vertices'): {path}")
    n = int(lines[1])
    verts: List[Tuple[float, float]] = []
    for i in range(n):
        x, y = map(float, lines[2 + i].split())
        verts.append((x, y))
    tri_hdr = 2 + n
    if lines[tri_hdr] != "# triangles":
        raise RuntimeError(f"Unexpected .tri format (missing '# triangles'): {path}")
    m = int(lines[tri_hdr + 1])
    tris: List[Tuple[int, int, int]] = []
    for i in range(m):
        a, b, c = map(int, lines[tri_hdr + 2 + i].split())
        tris.append((a, b, c))
    return Triangulation(verts=verts, tris=tris)


def write_svg(out_path: Path, title: str, tri: Triangulation, size_px: int = 900, margin_px: int = 30) -> None:
    pts = tri.verts
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    w = max(maxx - minx, 1e-9)
    h = max(maxy - miny, 1e-9)
    scale = (size_px - 2 * margin_px) / max(w, h)

    def tx(p: Tuple[float, float]) -> Tuple[float, float]:
        # SVG y axis goes down, so invert y.
        x, y = p
        sx = margin_px + (x - minx) * scale
        sy = margin_px + (maxy - y) * scale
        return sx, sy

    # Boundary path
    b = [tx(p) for p in pts]
    boundary_d = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in b) + " Z"

    # Unique undirected edges from triangles
    edges = set()
    for a, b_, c in tri.tris:
        for u, v in ((a, b_), (b_, c), (c, a)):
            if u == v:
                continue
            if u > v:
                u, v = v, u
            edges.add((u, v))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{size_px}" height="{size_px}" '
            f'viewBox="0 0 {size_px} {size_px}">\n'
        )
        f.write('<rect x="0" y="0" width="100%" height="100%" fill="white"/>\n')
        f.write(f'<text x="{margin_px}" y="{margin_px}" font-family="monospace" font-size="14">{title}</text>\n')

        # Triangle edges (light)
        f.write('<g stroke="#1f77b4" stroke-width="1" fill="none" opacity="0.85">\n')
        for u, v in sorted(edges):
            x1, y1 = tx(pts[u])
            x2, y2 = tx(pts[v])
            f.write(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}"/>\n')
        f.write("</g>\n")

        # Boundary on top
        f.write(f'<path d="{boundary_d}" stroke="#000000" stroke-width="2.5" fill="none"/>\n')

        f.write("</svg>\n")


def run_cli(exe: Path, poly_path: Path, tri_path: Path) -> None:
    proc = subprocess.run(
        [str(exe), "--input", str(poly_path), "--output", str(tri_path), "--algo", "chain_only"],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{exe.name} failed: rc={proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", type=Path, required=True, help="Path to reflex_cli binary")
    ap.add_argument("--outdir", type=Path, required=True, help="Directory to write SVGs")
    ap.add_argument("--n", type=int, default=200, help="Vertices for each example polygon")
    ap.add_argument("--seed", type=int, default=0, help="Seed for each family")
    args = ap.parse_args()

    exe = args.bin
    if not exe.exists():
        raise SystemExit(f"ERROR: binary not found: {exe}")

    tmp = args.outdir / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    families = [
        ("convex", convex_polygon),
        ("dent", dent_polygon),
        ("random", random_polygon),
    ]

    for name, gen in families:
        pts = gen(args.n, args.seed)
        poly_path = tmp / f"{name}.poly"
        tri_path = tmp / f"{name}.tri"
        write_poly(pts, poly_path)
        run_cli(exe, poly_path, tri_path)
        tri = read_tri(tri_path)
        out_svg = args.outdir / f"triangulation_{name}.svg"
        write_svg(out_svg, f"{name} (n={args.n}, seed={args.seed})", tri)


if __name__ == "__main__":
    main()

