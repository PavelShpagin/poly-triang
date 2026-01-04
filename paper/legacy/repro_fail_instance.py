#!/usr/bin/env python3
"""
Reproduce the first failing instance for the paper benchmark.

This is a debugging helper only:
- Generates polygons using the same generators as `paper/benchmark.py`
- Runs `build/bin/reflex_cli --algo chain`
- If the run fails, saves the polygon to `paper/_failing_<type>_<n>_seed<seed>.poly`
"""

from __future__ import annotations

import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, List, Tuple


ROOT = Path(__file__).resolve().parent.parent
EXE = ROOT / "build" / "bin" / "reflex_cli"
OUT_DIR = Path(__file__).resolve().parent

ROT = 0.123456789


def rotate_points(points: List[Tuple[float, float]], angle_rad: float) -> List[Tuple[float, float]]:
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in points]


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            f.write(f"{x:.17g} {y:.17g}\n")


def convex_polygon(n: int) -> List[Tuple[float, float]]:
    r = 100.0
    return [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]


def dent_polygon(n: int, dent_depth: float = 0.25) -> List[Tuple[float, float]]:
    outer = 100.0
    dent_i = max(0, n // 3)
    pts: List[Tuple[float, float]] = []
    for i in range(n):
        a = 2 * math.pi * i / n
        rr = outer * (dent_depth if i == dent_i else 1.0)
        pts.append((rr * math.cos(a), rr * math.sin(a)))
    return pts


def star_polygon(n: int) -> List[Tuple[float, float]]:
    n_pairs = max(3, n // 2)
    outer = 100.0
    inner = 30.0
    pts: List[Tuple[float, float]] = []
    for i in range(2 * n_pairs):
        angle = math.pi * i / n_pairs
        rr = outer if i % 2 == 0 else inner
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    if len(pts) > n:
        pts = pts[:n]
    elif len(pts) < n:
        base = pts[-1]
        for j in range(n - len(pts)):
            pts.append((base[0] + 1e-6 * (j + 1), base[1]))
    return pts


def random_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    rng = random.Random(seed + n)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    points: List[Tuple[float, float]] = []
    for angle in angles:
        rr = 100.0 * (0.4 + 0.6 * rng.random())
        points.append((rr * math.cos(angle), rr * math.sin(angle)))
    return points


def run_one(ptype: str, n: int, seed: int, gen: Callable[[int, int], List[Tuple[float, float]]]) -> bool:
    pts = rotate_points(gen(n, seed), ROT)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        poly_path = td / "poly.poly"
        out_path = td / "out.tri"
        write_poly(pts, poly_path)
        proc = subprocess.run(
            [str(EXE), "--input", str(poly_path), "--output", str(out_path), "--algo", "chain"],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            dest = OUT_DIR / f"_failing_{ptype}_{n}_seed{seed}.poly"
            write_poly(pts, dest)
            sys.stderr.write(f"FAIL: {ptype} n={n} seed={seed}\n")
            sys.stderr.write(proc.stderr.strip() + "\n")
            sys.stderr.write(f"Saved: {dest}\n")
            return False
    return True


def main() -> None:
    if not EXE.exists():
        raise SystemExit(f"Missing executable: {EXE}. Build first.")

    sizes = [100, 500, 1000, 2000, 5000, 10000]
    ppc = 5

    gens: dict[str, Callable[[int, int], List[Tuple[float, float]]]] = {
        "convex": lambda n, seed: convex_polygon(n),
        "dent": lambda n, seed: dent_polygon(n),
        "star": lambda n, seed: star_polygon(n),
        "random": lambda n, seed: random_polygon(n, seed),
    }

    for ptype, gen in gens.items():
        for n in sizes:
            for seed in range(ppc):
                ok = run_one(ptype, n, seed, gen)
                if not ok:
                    raise SystemExit(1)
    print("ALL OK")


if __name__ == "__main__":
    main()


