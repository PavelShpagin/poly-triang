#!/usr/bin/env python3
"""
Paper benchmark runner (4-way):
  - Ours (reflex_cli)              [C++]
  - Hertel-Mehlhorn (run_hertel.py) [Python baseline wrapper]
  - Garey et al. (run_garey.py)     [Python baseline wrapper]
  - Seidel / PolyTri (polytri_cli)  [C++]

Outputs:
  - paper/results_4way.csv

Notes:
  - This script is designed to be fast and robust for paper integration.
  - Sizes are intentionally capped to avoid extremely slow Python baselines.
  - No ear clipping baselines are included (per paper constraints).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple


@dataclass(frozen=True)
class Poly:
    path: Path
    polygon_name: str  # e.g. random_1000
    polygon_type: str  # e.g. random
    n: int
    r: int


TIME_RE = re.compile(r"time_ms=([0-9]+(?:\.[0-9]+)?)")
REFLEX_RE = re.compile(r"reflex_count=([0-9]+)")


def run_cmd(cmd: List[str], cwd: Path) -> str:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    out = proc.stdout
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{out}")
    return out


def parse_time_ms(output: str) -> float:
    m = TIME_RE.search(output)
    if not m:
        raise ValueError(f"Could not parse time_ms from output:\n{output}")
    return float(m.group(1))


def parse_reflex_count(output: str) -> int | None:
    m = REFLEX_RE.search(output)
    return int(m.group(1)) if m else None


def read_poly_points(poly_path: Path) -> List[Tuple[float, float]]:
    lines = poly_path.read_text(encoding="utf-8").strip().splitlines()
    n = int(lines[0].strip())
    pts: List[Tuple[float, float]] = []
    for i in range(1, 1 + n):
        x_str, y_str = lines[i].split()
        pts.append((float(x_str), float(y_str)))
    return pts


def signed_area(pts: List[Tuple[float, float]]) -> float:
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return a / 2.0


def count_reflex(pts: List[Tuple[float, float]]) -> int:
    n = len(pts)
    if n < 3:
        return 0
    ccw = signed_area(pts) > 0
    r = 0
    for i in range(n):
        ax, ay = pts[(i - 1) % n]
        bx, by = pts[i]
        cx, cy = pts[(i + 1) % n]
        cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        if ccw:
            if cross < -1e-10:
                r += 1
        else:
            if cross > 1e-10:
                r += 1
    return r


def discover_polys(poly_dir: Path, sizes: Iterable[int]) -> List[Poly]:
    polys: List[Poly] = []
    for p in sorted(poly_dir.glob("*.poly")):
        name = p.stem
        # pattern: type_size
        if "_" not in name:
            continue
        poly_type, n_str = name.rsplit("_", 1)
        try:
            n = int(n_str)
        except ValueError:
            continue
        if n not in set(sizes):
            continue
        pts = read_poly_points(p)
        r = count_reflex(pts)
        polys.append(Poly(path=p, polygon_name=name, polygon_type=poly_type, n=n, r=r))
    polys.sort(key=lambda x: (x.polygon_type, x.n))
    return polys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--sizes", nargs="+", type=int, default=[100, 500, 1000, 2000])
    ap.add_argument("--polygon-types", nargs="+", default=["convex", "random", "spiral", "star"])
    ap.add_argument("--runs", type=int, default=5)
    args = ap.parse_args()

    root = Path(args.project_root).resolve()
    poly_dir = root / "polygons" / "generated"
    out_csv = root / "paper" / "results_4way.csv"

    build_bin = root / "build" / "bin"
    reflex_cli = build_bin / "reflex_cli"
    polytri_cli = build_bin / "polytri_cli"

    if not reflex_cli.exists():
        raise RuntimeError(f"Missing {reflex_cli}. Build first.")
    if not polytri_cli.exists():
        raise RuntimeError(f"Missing {polytri_cli}. Build first.")

    # Ensure polygons exist (generate if needed)
    gen_script = root / "scripts" / "generate_polygons.py"
    if not poly_dir.exists() or not any(poly_dir.glob("*.poly")):
        print(f"[gen] generating polygons into {poly_dir}")
        run_cmd([sys.executable, str(gen_script), "--output", str(poly_dir), "--sizes", *map(str, args.sizes)], cwd=root)

    polys = discover_polys(poly_dir, args.sizes)
    polys = [p for p in polys if p.polygon_type in set(args.polygon_types)]
    if not polys:
        raise RuntimeError("No polygons found for requested sizes/types. Run generate_polygons.py first.")

    methods = [
        ("ours", lambda poly: run_cmd([str(reflex_cli), "--input", str(poly.path), "--output", str(root / "results" / f"reflex_{poly.polygon_name}.tri")], cwd=root)),
        ("hertel", lambda poly: run_cmd([sys.executable, str(root / "scripts" / "run_hertel.py"), "--input", str(poly.path), "--output", str(root / "results" / f"hertel_{poly.polygon_name}.tri")], cwd=root)),
        ("garey", lambda poly: run_cmd([sys.executable, str(root / "scripts" / "run_garey.py"), "--input", str(poly.path), "--output", str(root / "results" / f"garey_{poly.polygon_name}.tri")], cwd=root)),
        ("seidel_polytri", lambda poly: run_cmd([str(polytri_cli), "--input", str(poly.path), "--output", str(root / "results" / f"polytri_{poly.polygon_name}.tri")], cwd=root)),
    ]

    rows = []
    print("=== Paper 4-way Benchmark ===")
    print(f"root={root}")
    print(f"sizes={args.sizes}, polygon_types={args.polygon_types}, runs={args.runs}")
    print(f"polygons={len(polys)}")

    for poly in polys:
        print(f"\n[{poly.polygon_name}] n={poly.n} r={poly.r}")
        for alg, runner in methods:
            times = []
            for k in range(args.runs):
                out = runner(poly)
                t = parse_time_ms(out)
                times.append(t)
            times_sorted = sorted(times)
            median = times_sorted[len(times_sorted) // 2]
            rows.append({
                "algorithm": alg,
                "polygon": poly.polygon_name,
                "polygon_type": poly.polygon_type,
                "num_vertices": poly.n,
                "num_reflex": poly.r,
                "time_ms_median": median,
                "time_ms_runs": ";".join(f"{x:.6f}" for x in times),
            })
            print(f"  {alg:14s} median={median:10.4f} ms  runs={times}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["algorithm", "polygon", "polygon_type", "num_vertices", "num_reflex", "time_ms_median", "time_ms_runs"],
        )
        w.writeheader()
        w.writerows(rows)

    print(f"\nWrote {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


