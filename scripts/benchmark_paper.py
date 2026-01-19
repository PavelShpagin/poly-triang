#!/usr/bin/env python3
"""
Benchmark for O(n + r log r) polygon triangulation paper.

Compares our algorithm against classical baselines:
- Garey et al. (polypartition): O(n log n) monotone decomposition
- Hertel-Mehlhorn (polypartition): O(n log n) 
- Seidel: O(n log* n) randomized trapezoidal decomposition

Usage:
    python3 scripts/benchmark_paper.py [--sizes N1,N2,...] [--runs R]
"""

from __future__ import annotations
import argparse
import math
import os
import random
import statistics
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime

ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = ROOT / "build" / "bin"
RESULTS_DIR = ROOT / "results"

# Fixed rotation (radians) to avoid equal-y degeneracies in generated datasets.
ROT_ANGLE = 0.123456789

def log(msg: str) -> None:
    """Print timestamped log message."""
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}")


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    """Write polygon to file."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            # High precision to reduce accidental equal-y after rounding.
            f.write(f"{x:.17g} {y:.17g}\n")

def rotate_points(points: List[Tuple[float, float]], angle_rad: float) -> List[Tuple[float, float]]:
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in points]

def below_idx(points: List[Tuple[float, float]], a: int, b: int) -> bool:
    """Strict total order used for local-maxima counting (matches reflex_cli tie-break)."""
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
    """k = number of local maxima (equivalently, local minima) w.r.t sweep direction."""
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

def convex_polygon(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Convex polygon (k=1). Seeded affine image of a regular n-gon."""
    r = 100.0
    pts = [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]
    pts = rotate_points(pts, ROT_ANGLE)

    rng = random.Random(seed + 1337 * n)
    # Small seeded affine perturbation (keeps convexity).
    sx = 1.0 + 0.05 * (2.0 * rng.random() - 1.0)
    sy = 1.0 + 0.05 * (2.0 * rng.random() - 1.0)
    shx = 0.05 * (2.0 * rng.random() - 1.0)
    shy = 0.02 * (2.0 * rng.random() - 1.0)
    # Optional small extra rotation (still affine).
    ang = 0.15 * (2.0 * rng.random() - 1.0)
    pts = rotate_points(pts, ang)
    return [(sx * x + shx * y, shy * x + sy * y) for (x, y) in pts]

def dent_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    """Near-monotone polygon with a single dent (typically k=2)."""
    pts = convex_polygon(n, seed)
    if n < 5:
        return pts
    # Dent the current global maximum vertex (after the seeded affine), so we
    # reliably create exactly one extra extremum pair.
    imax = max(range(n), key=lambda i: (pts[i][1], -pts[i][0], -i))
    cx = sum(x for x, _ in pts) / n
    cy = sum(y for _, y in pts) / n

    rng = random.Random(seed + 4242 * n)
    depth = 0.55 + 0.25 * rng.random()  # [0.55, 0.80]

    def apply_depth(d: float) -> Tuple[float, float]:
        x, y = pts[imax]
        return (x + d * (cx - x), y + d * (cy - y))

    # Ensure dented vertex is strictly below its neighbors (in y-order), so it
    # stops being a local maximum and its neighbors become maxima -> k=2.
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
    """Random simple polygon via angular sweep - typically r ~ n/2."""
    rng = random.Random(seed + n)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    points = []
    for angle in angles:
        rr = 100.0 * (0.4 + 0.6 * rng.random())
        points.append((rr * math.cos(angle), rr * math.sin(angle)))
    return rotate_points(points, ROT_ANGLE)


def star_polygon(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Star polygon with r ~ n/2 reflex vertices."""
    rng = random.Random(seed + 9001 * n)
    n_pairs = max(3, n // 2)
    # Seeded radii + phase (mild jitter).
    outer = 100.0 * (0.92 + 0.16 * rng.random())
    inner = 30.0 * (0.90 + 0.20 * rng.random())
    phase = 2 * math.pi * rng.random()
    pts = []
    for i in range(2 * n_pairs):
        angle = phase + math.pi * i / n_pairs
        rr = outer if i % 2 == 0 else inner
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    while len(pts) > n:
        pts.pop()
    while len(pts) < n:
        pts.append((pts[-1][0] + 1e-7 * (len(pts) - n + 1), pts[-1][1]))
    return rotate_points(pts, ROT_ANGLE)


def spiral_polygon(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Spiral polygon with r ~ 1 reflex vertices (near-convex)."""
    turns = 3.0
    start, end = 20.0, 100.0
    pts = []
    for i in range(n):
        t = i / max(1, n - 1)
        angle = 2 * math.pi * turns * t
        rr = start + (end - start) * t
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    return pts


@dataclass
class RunResult:
    time_ms: float
    triangles: int
    reflex_count: Optional[int] = None


def parse_kv_output(stdout: str) -> Dict[str, str]:
    """Parse key=value output from CLI."""
    line = stdout.strip().splitlines()[-1].strip() if stdout.strip() else ""
    parts = line.split(",")
    out: Dict[str, str] = {"algorithm": parts[0] if parts else ""}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def run_cli(exe: Path, poly_path: Path, out_path: Path, timeout: int = 300) -> Optional[RunResult]:
    """Run triangulation CLI and parse results. Handles crashes gracefully."""
    if not exe.exists():
        return None
    try:
        cmd = [str(exe), "--input", str(poly_path), "--output", str(out_path)]
        # Be explicit about implementation choice for our binary (no hidden fallbacks).
        if exe.name == "reflex_cli":
            cmd += ["--algo", "chain"]
        proc = subprocess.run(
            cmd,
            capture_output=True, text=True, timeout=timeout
        )
        if proc.returncode != 0:
            log(f"ERROR: {exe.name} failed (code={proc.returncode})")
            if proc.stdout.strip():
                print(proc.stdout, end="" if proc.stdout.endswith("\n") else "\n")
            if proc.stderr.strip():
                print(proc.stderr, end="" if proc.stderr.endswith("\n") else "\n", file=sys.stderr)
            return None
        kv = parse_kv_output(proc.stdout)
        time_ms = float(kv.get("time_ms", "0"))
        triangles = int(kv.get("triangles", "0"))
        reflex = int(kv["reflex_count"]) if "reflex_count" in kv else None
        return RunResult(time_ms=time_ms, triangles=triangles, reflex_count=reflex)
    except subprocess.TimeoutExpired:
        log(f"ERROR: {exe.name} timed out after {timeout}s")
        return None
    except Exception:
        log(f"ERROR: {exe.name} crashed (exception)")
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark polygon triangulation")
    parser.add_argument("--sizes", default="1000,10000,100000",
                       help="Comma-separated polygon sizes (default: 1000,10000,100000)")
    parser.add_argument("--runs", type=int, default=10,
                       help="Number of runs per configuration (default: 10)")
    parser.add_argument("--types", default="convex,random,star,spiral",
                       help="Comma-separated polygon families to benchmark (default: convex,random,star,spiral)")
    parser.add_argument("--out-csv", default=str(RESULTS_DIR / "benchmark_results.csv"),
                       help="Output CSV path for aggregated results (default: results/benchmark_results.csv)")
    parser.add_argument("--out-raw-csv", default=str(RESULTS_DIR / "benchmark_results_raw.csv"),
                       help="Output CSV path for per-run raw results (default: results/benchmark_results_raw.csv)")
    parser.add_argument("--timeout", type=int, default=300,
                       help="Per-run timeout in seconds (default: 300)")
    args = parser.parse_args()
    
    sizes = [int(s.strip()) for s in args.sizes.split(",")]
    runs_per_config = args.runs
    
    # Algorithms to compare
    executables = {
        "ours": BIN_DIR / "reflex_cli",
        "garey": BIN_DIR / "polypartition_mono_cli",
        "hertel": BIN_DIR / "polypartition_hm_cli",
        "seidel": BIN_DIR / "seidel_cli",
    }
    
    # Check which executables exist
    available = {k: p for k, p in executables.items() if p.exists()}
    missing = [k for k, p in executables.items() if not p.exists()]
    
    if "ours" not in available:
        log(f"ERROR: reflex_cli not found!")
        log(f"Build with: cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j")
        sys.exit(1)
    
    if missing:
        log(f"Note: Some baselines not available: {missing}")
    
    polygon_types = {
        "convex": convex_polygon,
        "dent": dent_polygon,
        "random": random_polygon,
        "star": star_polygon,
        "spiral": spiral_polygon,
    }

    selected_types = [t.strip() for t in args.types.split(",") if t.strip()]
    unknown = sorted(set(selected_types) - set(polygon_types.keys()))
    if unknown:
        log(f"ERROR: Unknown polygon types: {unknown} (known: {list(polygon_types.keys())})")
        sys.exit(2)
    polygon_types = {k: v for k, v in polygon_types.items() if k in set(selected_types)}
    
    log(f"Starting benchmark: sizes={sizes}, runs={runs_per_config}")
    log(f"Algorithms: {list(available.keys())}")
    
    results = {}
    RESULTS_DIR.mkdir(exist_ok=True)
    out_csv_path = Path(args.out_csv)
    out_raw_csv_path = Path(args.out_raw_csv)
    out_csv_path.parent.mkdir(parents=True, exist_ok=True)
    out_raw_csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Aggregated CSV (kept for backward compatibility with earlier ad-hoc scripts).
    with open(out_csv_path, "w", encoding="utf-8") as csv_file:
        header = "polygon_type,n,r,ours_ms,garey_ms,hertel_ms,seidel_ms\n"
        csv_file.write(header)

        # Raw per-run CSV (used to regenerate CGAT LaTeX tables with mean±stdev).
        with open(out_raw_csv_path, "w", encoding="utf-8") as raw_csv:
            raw_csv.write("polygon_type,n,seed,k_count,algorithm,time_ms,triangles,reflex_count\n")
        
            with tempfile.TemporaryDirectory() as td:
                td = Path(td)
                poly_path = td / "poly.poly"
                out_path = td / "out.tri"
                
                for ptype, gen in polygon_types.items():
                    log(f"\n=== {ptype.upper()} polygons ===")
                    
                    for n in sizes:
                        times: Dict[str, List[float]] = {k: [] for k in available.keys()}
                        reflex_vals: List[int] = []
                        k_vals: List[int] = []
                        
                        for seed in range(runs_per_config):
                            pts = gen(n, seed)
                            k_count = count_local_maxima_k(pts)
                            k_vals.append(k_count)
                            write_poly(pts, poly_path)
                            
                            # Run each algorithm
                            for alg, exe in available.items():
                                result = run_cli(exe, poly_path, out_path, timeout=args.timeout)
                                if result and result.triangles == n - 2:
                                    times[alg].append(result.time_ms)
                                    if alg == "ours" and result.reflex_count is not None:
                                        reflex_vals.append(result.reflex_count)
                                    raw_csv.write(
                                        f"{ptype},{n},{seed},{k_count},{alg},{result.time_ms},{result.triangles},{result.reflex_count or ''}\n"
                                    )
                                else:
                                    raw_csv.write(f"{ptype},{n},{seed},{k_count},{alg},,,\n")
                        
                        r_avg = statistics.mean(reflex_vals) if reflex_vals else 0.0
                        means: Dict[str, Optional[float]] = {}
                        for alg, vals in times.items():
                            means[alg] = statistics.mean(vals) if vals else None

                        # Calculate speedups vs baselines
                        ours_mean = means.get("ours")
                        speedups = {}
                        for alg in ["garey", "hertel", "seidel"]:
                            if ours_mean and means.get(alg) and ours_mean > 0:
                                speedups[alg] = means[alg] / ours_mean

                        # Log progress
                        ours_str = f"{ours_mean:.3f}" if ours_mean else "--"
                        garey_str = f"{means.get('garey', 0):.3f}" if means.get('garey') else "--"
                        hertel_str = f"{means.get('hertel', 0):.3f}" if means.get('hertel') else "--"
                        seidel_str = f"{means.get('seidel', 0):.3f}" if means.get('seidel') else "--"

                        best_speedup = max(speedups.values()) if speedups else 0
                        speedup_str = f"{best_speedup:.1f}x" if best_speedup > 0 else "--"

                        log(f"  n={n:>7,} r={int(r_avg):>6,} ({r_avg/n*100:4.0f}%): ours={ours_str:>8}ms garey={garey_str:>8}ms seidel={seidel_str:>8}ms [{speedup_str}]")

                        # Save to results dict
                        results[(ptype, n)] = {
                            "type": ptype, "n": n, "r_avg": r_avg,
                            "ours_mean": ours_mean,
                            "garey_mean": means.get("garey"),
                            "hertel_mean": means.get("hertel"),
                            "seidel_mean": means.get("seidel"),
                            "speedups": speedups
                        }

                        # Write to CSV
                        csv_file.write(f"{ptype},{n},{int(r_avg)},{ours_mean or ''},{means.get('garey') or ''},{means.get('hertel') or ''},{means.get('seidel') or ''}\n")
    
    log(f"\nResults saved to {out_csv_path}")
    log(f"Raw per-run results saved to {out_raw_csv_path}")
    
    # Generate summary
    log("\n" + "=" * 70)
    log("SUMMARY: Our O(n + r log r) vs O(n log n) baselines")
    log("=" * 70)
    log(f"{'Type':<10} {'n':>10} {'r':>8} {'Ours (ms)':>12} {'Garey':>12} {'Seidel':>12}")
    log("-" * 70)
    
    for ptype in polygon_types.keys():
        for n in sizes:
            rec = results.get((ptype, n))
            if rec:
                ours = f"{rec['ours_mean']:.3f}" if rec['ours_mean'] else "--"
                garey = f"{rec['garey_mean']:.3f}" if rec.get('garey_mean') else "--"
                seidel = f"{rec['seidel_mean']:.3f}" if rec.get('seidel_mean') else "--"
                log(f"{ptype:<10} {n:>10,} {int(rec['r_avg']):>8,} {ours:>12} {garey:>12} {seidel:>12}")
    
    log("\nKey insight: For r << n, our O(n + r log r) beats O(n log n).")
    log("Benchmark complete!")


if __name__ == "__main__":
    main()

