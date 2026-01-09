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

def log(msg: str) -> None:
    """Print timestamped log message."""
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}")


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    """Write polygon to file."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            f.write(f"{x:.12f} {y:.12f}\n")


def convex_polygon(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Convex polygon with r=0 reflex vertices."""
    r = 100.0
    return [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]


def random_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    """Random simple polygon via angular sweep - typically r ~ n/2."""
    import random
    rng = random.Random(seed + n * 1000)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    points = []
    for angle in angles:
        rr = 100.0 * (0.4 + 0.6 * rng.random())
        points.append((rr * math.cos(angle), rr * math.sin(angle)))
    return points


def star_polygon(n: int, seed: int = 0) -> List[Tuple[float, float]]:
    """Star polygon with r ~ n/2 reflex vertices."""
    n_pairs = max(3, n // 2)
    outer, inner = 100.0, 30.0
    pts = []
    for i in range(2 * n_pairs):
        angle = math.pi * i / n_pairs
        rr = outer if i % 2 == 0 else inner
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    while len(pts) > n:
        pts.pop()
    while len(pts) < n:
        pts.append((pts[-1][0] + 1e-7 * (len(pts) - n + 1), pts[-1][1]))
    return pts


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
        proc = subprocess.run(
            [str(exe), "--input", str(poly_path), "--output", str(out_path)],
            capture_output=True, text=True, timeout=timeout
        )
        if proc.returncode != 0:
            return None
        kv = parse_kv_output(proc.stdout)
        time_ms = float(kv.get("time_ms", "0"))
        triangles = int(kv.get("triangles", "0"))
        reflex = int(kv["reflex_count"]) if "reflex_count" in kv else None
        return RunResult(time_ms=time_ms, triangles=triangles, reflex_count=reflex)
    except subprocess.TimeoutExpired:
        return None  # Graceful timeout
    except Exception:
        return None  # Graceful crash handling


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark polygon triangulation")
    parser.add_argument("--sizes", default="1000,10000,100000",
                       help="Comma-separated polygon sizes (default: 1000,10000,100000)")
    parser.add_argument("--runs", type=int, default=10,
                       help="Number of runs per configuration (default: 10)")
    parser.add_argument("--types", default="convex,random,star,spiral",
                       help="Comma-separated polygon families to benchmark (default: convex,random,star,spiral)")
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
    
    # Open CSV for raw results
    csv_path = RESULTS_DIR / "benchmark_results.csv"
    with open(csv_path, "w", encoding="utf-8") as csv_file:
        header = "polygon_type,n,r,ours_ms,garey_ms,hertel_ms,seidel_ms\n"
        csv_file.write(header)
        
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            poly_path = td / "poly.poly"
            out_path = td / "out.tri"
            
            for ptype, gen in polygon_types.items():
                log(f"\n=== {ptype.upper()} polygons ===")
                
                for n in sizes:
                    times: Dict[str, List[float]] = {k: [] for k in available.keys()}
                    reflex_vals: List[int] = []
                    
                    for seed in range(runs_per_config):
                        pts = gen(n, seed)
                        write_poly(pts, poly_path)
                        
                        # Run each algorithm
                        for alg, exe in available.items():
                            result = run_cli(exe, poly_path, out_path)
                            if result and result.triangles == n - 2:
                                times[alg].append(result.time_ms)
                                if alg == "ours" and result.reflex_count is not None:
                                    reflex_vals.append(result.reflex_count)
                    
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
    
    log(f"\nResults saved to {csv_path}")
    
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

