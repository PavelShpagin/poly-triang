#!/usr/bin/env python3
"""
Self-contained benchmark harness for the CGAT artifact.

Writes:
  - generated/benchmark_results_raw.csv  (per-run)
  - generated/benchmark_results.csv      (aggregated means, for convenience)
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parent.parent
ROT_ANGLE = 0.123456789


def rotate_points(points: List[Tuple[float, float]], angle_rad: float) -> List[Tuple[float, float]]:
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in points]


def below_idx(points: List[Tuple[float, float]], a: int, b: int) -> bool:
    ay = points[a][1]
    by = points[b][1]
    if ay < by:
        return True
    if ay > by:
        return False
    ax = points[a][0]
    bx = points[b][0]
    if ax > bx:
        return True
    if ax < bx:
        return False
    return a > b


def count_local_maxima_k(points: List[Tuple[float, float]]) -> int:
    n = len(points)
    if n < 3:
        return 0
    k = 0
    for i in range(n):
        p = (i - 1 + n) % n
        nx = (i + 1) % n
        if below_idx(points, p, i) and below_idx(points, nx, i):
            k += 1
    return k


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


def star_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    rng = random.Random(seed + 9001 * n)
    n_pairs = max(3, n // 2)
    outer = 100.0 * (0.92 + 0.16 * rng.random())
    inner = 30.0 * (0.90 + 0.20 * rng.random())
    phase = 2 * math.pi * rng.random()
    pts: List[Tuple[float, float]] = []
    for i in range(2 * n_pairs):
        angle = phase + math.pi * i / n_pairs
        rr = outer if i % 2 == 0 else inner
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    while len(pts) > n:
        pts.pop()
    while len(pts) < n:
        pts.append((pts[-1][0] + 1e-7 * (len(pts) - n + 1), pts[-1][1]))
    return rotate_points(pts, ROT_ANGLE)


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            f.write(f"{x:.17g} {y:.17g}\n")


def parse_kv_output(stdout: str) -> Dict[str, str]:
    line = stdout.strip().splitlines()[-1].strip() if stdout.strip() else ""
    parts = line.split(",")
    out: Dict[str, str] = {"algorithm": parts[0] if parts else ""}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k.strip()] = v.strip()
    return out


@dataclass
class RunResult:
    time_ms: float
    triangles: int
    reflex_count: Optional[int] = None


def run_cli(exe: Path, poly_path: Path, out_path: Path, timeout: int) -> Optional[RunResult]:
    if not exe.exists():
        return None
    try:
        proc = subprocess.run(
            [str(exe), "--input", str(poly_path), "--output", str(out_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if proc.returncode != 0:
            return None
        kv = parse_kv_output(proc.stdout)
        time_ms = float(kv.get("time_ms", "0"))
        triangles = int(kv.get("triangles", "0"))
        reflex = int(kv["reflex_count"]) if "reflex_count" in kv else None
        return RunResult(time_ms=time_ms, triangles=triangles, reflex_count=reflex)
    except subprocess.TimeoutExpired:
        return None
    except Exception:
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="100,500,1000,2000,5000,10000")
    ap.add_argument("--runs", type=int, default=5, help="Seeds 0..runs-1 per config")
    ap.add_argument("--timeout", type=int, default=10)
    ap.add_argument("--types", default="convex,dent,random,star")
    ap.add_argument("--bin-dir", type=Path, default=ROOT / "bin")
    ap.add_argument("--out-raw", type=Path, default=ROOT / "generated" / "benchmark_results_raw.csv")
    ap.add_argument("--out-csv", type=Path, default=ROOT / "generated" / "benchmark_results.csv")
    args = ap.parse_args()

    sizes = [int(s.strip()) for s in args.sizes.split(",") if s.strip()]
    runs = args.runs
    timeout = args.timeout

    generators = {
        "convex": convex_polygon,
        "dent": dent_polygon,
        "random": random_polygon,
        "star": star_polygon,
    }
    selected = [t.strip() for t in args.types.split(",") if t.strip()]
    unknown = sorted(set(selected) - set(generators.keys()))
    if unknown:
        print(f"ERROR: unknown types: {unknown}", file=sys.stderr)
        sys.exit(2)
    generators = {k: v for k, v in generators.items() if k in set(selected)}

    exes = {
        "ours": args.bin_dir / "reflex_cli",
        "garey": args.bin_dir / "polypartition_mono_cli",
        "hertel": args.bin_dir / "polypartition_hm_cli",
        "seidel": args.bin_dir / "seidel_cli",
    }

    args.out_raw.parent.mkdir(parents=True, exist_ok=True)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        poly_path = td / "poly.poly"
        out_path = td / "out.tri"

        # Raw CSV
        with open(args.out_raw, "w", encoding="utf-8", newline="") as fraw:
            wraw = csv.writer(fraw)
            wraw.writerow(["polygon_type", "n", "seed", "k_count", "algorithm", "time_ms", "triangles", "reflex_count"])

            # Aggregate data
            agg_times: Dict[Tuple[str, int, str], List[float]] = {}
            agg_reflex: Dict[Tuple[str, int], List[int]] = {}

            for ptype, gen in generators.items():
                for n in sizes:
                    for seed in range(runs):
                        pts = gen(n, seed)
                        k_count = count_local_maxima_k(pts)
                        write_poly(pts, poly_path)

                        for alg, exe in exes.items():
                            res = run_cli(exe, poly_path, out_path, timeout=timeout)
                            if res and res.triangles == n - 2:
                                wraw.writerow([ptype, n, seed, k_count, alg, f"{res.time_ms:.6f}", res.triangles, res.reflex_count or ""])
                                agg_times.setdefault((ptype, n, alg), []).append(res.time_ms)
                                if alg == "ours" and res.reflex_count is not None:
                                    agg_reflex.setdefault((ptype, n), []).append(res.reflex_count)
                            else:
                                wraw.writerow([ptype, n, seed, k_count, alg, "", "", ""])

        # Aggregated CSV (means only, for convenience)
        with open(args.out_csv, "w", encoding="utf-8", newline="") as fout:
            w = csv.writer(fout)
            w.writerow(["polygon_type", "n", "r", "ours_ms", "garey_ms", "hertel_ms", "seidel_ms"])
            for ptype in generators.keys():
                for n in sizes:
                    r_vals = agg_reflex.get((ptype, n), [])
                    r_mean = int(round(statistics.mean(r_vals))) if r_vals else 0
                    row = [ptype, n, r_mean]
                    for alg in ["ours", "garey", "hertel", "seidel"]:
                        vals = agg_times.get((ptype, n, alg), [])
                        row.append(f"{statistics.mean(vals):.6f}" if vals else "")
                    w.writerow(row)


if __name__ == "__main__":
    main()

