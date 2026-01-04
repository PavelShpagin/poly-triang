#!/usr/bin/env python3
"""
Benchmark: Ours (reflex_cli) vs Seidel (polytri_cli).

- Uses CLI-reported time_ms (no subprocess timing / I/O overhead)
- Generates multiple polygons per (type, n) and reports mean ± std.
"""

import math
import random
import statistics
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
BIN_DIR = ROOT / "build" / "bin"


def gen_random_star(n: int, seed: int):
    random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = random.uniform(30, 100)
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_convex(n: int, seed: int):
    # Deterministic convex polygon (seed ignored to keep r=0 always).
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100.0
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_star(n: int, seed: int):
    # Deterministic star polygon (seed ignored).
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100.0 if i % 2 == 0 else 40.0
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_spiral(n: int, seed: int):
    # Deterministic spiral polygon (seed ignored).
    pts = []
    for i in range(n):
        t = i / n
        angle = 4 * math.pi * t
        r = 20 + 80 * t
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def write_polygon(pts, path: Path):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")


def run_cli(cli: Path, input_path: Path, output_path: Path):
    try:
        proc = subprocess.run(
            [str(cli), "-i", str(input_path), "-o", str(output_path)],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout


def parse_time_ms(stdout: str):
    for line in stdout.splitlines():
        if "time_ms=" in line:
            return float(line.split("time_ms=")[1].split(",")[0].split()[0])
    return None


def parse_reflex_count(stdout: str):
    for line in stdout.splitlines():
        if "reflex_count=" in line:
            return int(line.split("reflex_count=")[1].split(",")[0].split()[0])
    return None


def main():
    reflex_cli = BIN_DIR / "reflex_cli"
    polytri_cli = BIN_DIR / "polytri_cli"
    if not reflex_cli.exists():
        raise SystemExit("Missing build/bin/reflex_cli")
    if not polytri_cli.exists():
        raise SystemExit("Missing build/bin/polytri_cli")

    polygon_types = {
        "convex": gen_convex,
        "spiral": gen_spiral,
        "random": gen_random_star,
        "star": gen_star,
    }
    sizes = [1_000, 5_000, 10_000, 20_000, 50_000, 100_000]
    polygons_per_config = 30

    print("=" * 96)
    print("BENCHMARK: Ours (reflex_cli) vs Seidel (polytri_cli)")
    print(f"polygons_per_config={polygons_per_config}")
    print("=" * 96)

    results = {}

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        in_path = td / "poly.poly"
        out_path = td / "out.tri"

        for ptype, gen in polygon_types.items():
            print(f"\n{ptype.upper()}:")
            for n in sizes:
                ours = []
                seidel = []
                rvals = []

                for seed in range(polygons_per_config):
                    pts = gen(n, seed)
                    write_polygon(pts, in_path)

                    out_reflex = run_cli(reflex_cli, in_path, out_path)
                    out_polytri = run_cli(polytri_cli, in_path, out_path)
                    if out_reflex is None or out_polytri is None:
                        continue

                    t_reflex = parse_time_ms(out_reflex)
                    t_polytri = parse_time_ms(out_polytri)
                    r = parse_reflex_count(out_reflex)
                    if t_reflex is None or t_polytri is None or r is None:
                        continue

                    ours.append(t_reflex)
                    seidel.append(t_polytri)
                    rvals.append(r)

                if not ours:
                    continue

                ours_mean = statistics.mean(ours)
                ours_std = statistics.stdev(ours) if len(ours) > 1 else 0.0
                seid_mean = statistics.mean(seidel)
                seid_std = statistics.stdev(seidel) if len(seidel) > 1 else 0.0
                r_mean = statistics.mean(rvals)
                speedup = seid_mean / ours_mean if ours_mean > 0 else float("inf")

                results[(ptype, n)] = {
                    "r_mean": r_mean,
                    "r_ratio": r_mean / n,
                    "ours_mean": ours_mean,
                    "ours_std": ours_std,
                    "seid_mean": seid_mean,
                    "seid_std": seid_std,
                    "speedup": speedup,
                    "samples": len(ours),
                }

                print(
                    f"  n={n:>7}  r_avg={int(r_mean):>6} ({100*r_mean/n:>5.1f}%)  "
                    f"Ours={ours_mean:>8.2f}±{ours_std:>6.2f}ms  "
                    f"Seidel={seid_mean:>8.2f}±{seid_std:>6.2f}ms  "
                    f"Speedup={speedup:>6.2f}x"
                )

    print("\n" + "=" * 96)
    print("LATEX (main table, n=100000)")
    print("=" * 96)
    print(r"\\begin{table}[t]")
    print(r"\\centering")
    print(r"\\caption{Running times (ms) at $n = 100{,}000$ vertices, averaged over 30 polygons.}")
    print(r"\\label{tab:results}")
    print(r"\\smallskip")
    print(r"\\begin{tabular}{lrrrrr}")
    print(r"\\toprule")
    print(r"\\textbf{Type} & $r$ (avg) & $r/n$ & \\textbf{Ours} & \\textbf{Seidel} & \\textbf{Speedup} \\\\")
    print(r"\\midrule")
    for ptype in ["convex", "spiral", "random", "star"]:
        key = (ptype, 100_000)
        if key not in results:
            continue
        r = results[key]
        print(
            f"{ptype.title()} & {int(r['r_mean']):,} & {100*r['r_ratio']:.0f}\\% & "
            f"${r['ours_mean']:.1f} \\pm {r['ours_std']:.1f}$ & "
            f"${r['seid_mean']:.1f} \\pm {r['seid_std']:.1f}$ & "
            f"{r['speedup']:.1f}$\\times$ \\\\"
        )
    print(r"\\bottomrule")
    print(r"\\end{tabular}")
    print(r"\\end{table}")


if __name__ == "__main__":
    main()


