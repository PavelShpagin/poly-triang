#!/usr/bin/env python3
"""
Paper benchmark runner (final).

Compares:
- ours:            build/bin/reflex_cli
- garey baseline:  build/bin/polypartition_mono_cli   (monotone partition + triangulation)
- hertel baseline: build/bin/polypartition_hm_cli     (HM convex partition + triangulate pieces)
- seidel baseline: build/bin/seidel_cli               (palmerc/Seidel C implementation)

Outputs:
- paper/results.txt                 (human-readable)
- paper/generated/benchmark_table.tex (random polygons table)
- paper/generated/benchmark_full.tex  (all families table)

Notes:
- We use CLI-reported time_ms for the C++ baselines that print it.
- For seidel_cli we use its printed time_ms.
- Polygon generator matches scripts/generate_polygons.py (with multiple seeds for random).
"""

from __future__ import annotations

import argparse
import math
import statistics
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = ROOT / "build" / "bin"
OUT_TXT = Path(__file__).resolve().parent / "results.txt"
OUT_DIR = Path(__file__).resolve().parent / "generated"


def write_poly(points: List[Tuple[float, float]], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            # Use high precision to avoid accidental degeneracies after rounding.
            f.write(f"{x:.17g} {y:.17g}\n")


def rotate_points(points: List[Tuple[float, float]], angle_rad: float) -> List[Tuple[float, float]]:
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return [(ca * x - sa * y, sa * x + ca * y) for (x, y) in points]


def convex_polygon(n: int) -> List[Tuple[float, float]]:
    r = 100.0
    return [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]


def random_polygon(n: int, seed: int) -> List[Tuple[float, float]]:
    # Same model as scripts/generate_polygons.py
    import random

    rng = random.Random(seed + n)
    angles = sorted(rng.random() * 2 * math.pi for _ in range(n))
    points = []
    for angle in angles:
        rr = 100.0 * (0.4 + 0.6 * rng.random())
        points.append((rr * math.cos(angle), rr * math.sin(angle)))
    return points


def star_polygon(n: int) -> List[Tuple[float, float]]:
    n_pairs = max(3, n // 2)
    outer = 100.0
    inner = 30.0
    pts = []
    for i in range(2 * n_pairs):
        angle = math.pi * i / n_pairs
        rr = outer if i % 2 == 0 else inner
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    # Ensure exactly n vertices
    if len(pts) > n:
        pts = pts[:n]
    elif len(pts) < n:
        # pad by repeating last point slightly perturbed in angle (keeps simple enough for our usage)
        k = n - len(pts)
        base = pts[-1]
        for j in range(k):
            pts.append((base[0] + 1e-6 * (j + 1), base[1]))
    return pts


def spiral_polygon(n: int) -> List[Tuple[float, float]]:
    turns = 3.0
    start = 20.0
    end = 100.0
    pts = []
    for i in range(n):
        t = i / (n - 1)
        angle = 2 * math.pi * turns * t
        rr = start + (end - start) * t
        pts.append((rr * math.cos(angle), rr * math.sin(angle)))
    return pts


def dent_polygon(n: int, dent_depth: float = 0.25) -> List[Tuple[float, float]]:
    # Simple star-shaped polygon with (typically) a very small number of reflex vertices:
    # a nearly circular polygon with a single inward "dent".
    pts = []
    outer = 100.0
    dent_i = max(0, n // 3)
    for i in range(n):
        a = 2 * math.pi * i / n
        rr = outer * (dent_depth if i == dent_i else 1.0)
        pts.append((rr * math.cos(a), rr * math.sin(a)))
    return pts


@dataclass
class RunResult:
    time_ms: float
    triangles: int
    reflex_count: int | None = None


def parse_kv_line(stdout: str) -> Dict[str, str]:
    # expects: name,k=v,k=v,...
    line = stdout.strip().splitlines()[-1].strip()
    parts = line.split(",")
    out: Dict[str, str] = {}
    out["algorithm"] = parts[0]
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def run_cli(exe: Path, poly_path: Path, out_path: Path, timeout_s: float) -> RunResult:
    proc = subprocess.run(
        [str(exe), "--input", str(poly_path), "--output", str(out_path)],
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"{exe.name} failed")
    kv = parse_kv_line(proc.stdout)
    time_ms = float(kv["time_ms"])
    triangles = int(kv["triangles"])
    reflex = int(kv["reflex_count"]) if "reflex_count" in kv else None
    return RunResult(time_ms=time_ms, triangles=triangles, reflex_count=reflex)

def safe_run_cli(exe: Path, poly_path: Path, out_path: Path, timeout_s: float) -> Optional[RunResult]:
    try:
        return run_cli(exe, poly_path, out_path, timeout_s=timeout_s)
    except Exception:
        return None

def fmt_pm(mean: Optional[float], std: Optional[float]) -> str:
    if mean is None:
        return "--"
    if std is None:
        return f"{mean:.2f}"
    return f"{mean:.2f} $\\pm$ {std:.2f}"

def write_tex_tables(summary: Dict[Tuple[str, int], dict], sizes: List[int], polygon_types: List[str]) -> None:
    OUT_DIR.mkdir(exist_ok=True)

    # Main table: random polygons
    lines: List[str] = []
    lines.append(r"\begin{table}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{Running time (ms) on random polygons (mean $\pm$ stdev over instances).}")
    lines.append(r"\label{tab:benchmark}")
    lines.append(r"\begin{tabular}{rrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"$n$ & $r$ & \textbf{Ours} & Seidel & Garey & Hertel--Mehlhorn \\")
    lines.append(r"\midrule")
    for n in sizes:
        rec = summary.get(("random", n))
        if not rec:
            continue
        r_avg = int(rec["r_avg"])
        cols = [
            ("ours", rec.get("ours_mean"), rec.get("ours_std")),
            ("seidel", rec.get("seidel_mean"), rec.get("seidel_std")),
            ("garey", rec.get("garey_mean"), rec.get("garey_std")),
            ("hertel", rec.get("hertel_mean"), rec.get("hertel_std")),
        ]
        # Bold the fastest available mean in row
        best = min((m for _, m, _ in cols if m is not None), default=None)
        vals = []
        for _, m, s in cols:
            cell = fmt_pm(m, s)
            if best is not None and m is not None and abs(m - best) < 1e-12:
                cell = r"\textbf{" + cell + "}"
            vals.append(cell)
        lines.append(f"{n:,} & {r_avg} & " + " & ".join(vals) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    (OUT_DIR / "benchmark_table.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Full table: all families
    lines = []
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{Running time (ms) across polygon families (mean over instances).}")
    lines.append(r"\label{tab:benchmark-full}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{llrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"Type & $n$ & $r$ & \textbf{Ours} & Seidel & Garey & Hertel--Mehlhorn \\")
    lines.append(r"\midrule")
    for ptype in polygon_types:
        first = True
        for n in sizes:
            rec = summary.get((ptype, n))
            if not rec:
                continue
            r_avg = int(rec["r_avg"])
            cols = [
                ("ours", rec.get("ours_mean"), rec.get("ours_std")),
                ("seidel", rec.get("seidel_mean"), rec.get("seidel_std")),
                ("garey", rec.get("garey_mean"), rec.get("garey_std")),
                ("hertel", rec.get("hertel_mean"), rec.get("hertel_std")),
            ]
            best = min((m for _, m, _ in cols if m is not None), default=None)
            vals = []
            for _, m, s in cols:
                cell = fmt_pm(m, s)
                if best is not None and m is not None and abs(m - best) < 1e-12:
                    cell = r"\textbf{" + cell + "}"
                vals.append(cell)
            type_col = ptype.capitalize() if first else ""
            first = False
            lines.append(f"{type_col} & {n:,} & {r_avg} & " + " & ".join(vals) + r" \\")
        lines.append(r"\midrule")
    if lines[-1] == r"\midrule":
        lines[-1] = r"\bottomrule"
    else:
        lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table*}")
    (OUT_DIR / "benchmark_full.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Run paper benchmark with robust timeouts.")
    ap.add_argument("--sizes", nargs="+", type=int, default=[100, 500, 1000, 2000, 5000, 10000])
    ap.add_argument("--polygons-per-config", type=int, default=5,
                    help="Number of random instances per (type,n). Keep small for quick runs.")
    ap.add_argument("--timeout-ours", type=float, default=5.0)
    ap.add_argument("--timeout-garey", type=float, default=5.0)
    ap.add_argument("--timeout-hertel", type=float, default=5.0)
    ap.add_argument("--timeout-seidel", type=float, default=5.0)
    ap.add_argument("--skip-seidel", action="store_true", help="Skip seidel baseline entirely.")
    args = ap.parse_args()

    # Verify executables exist
    required = {
        "ours": BIN_DIR / "reflex_cli",
        "garey": BIN_DIR / "polypartition_mono_cli",
        "hertel": BIN_DIR / "polypartition_hm_cli",
        "seidel": BIN_DIR / "seidel_cli",
    }
    if args.skip_seidel:
        required.pop("seidel", None)
    missing = [k for k, p in required.items() if not p.exists()]
    if missing:
        raise SystemExit(f"Missing executables: {missing}. Build first.")

    # NOTE: some baselines may reject spiral polygons if degeneracies occur;
    # benchmark.py will record failures as missing in the TeX tables.
    # Deterministic rotation: helps keep inputs in general position for sweep algorithms.
    rot = 0.123456789
    polygon_types = {
        "convex": lambda n, seed: rotate_points(convex_polygon(n), rot),
        "dent": lambda n, seed: rotate_points(dent_polygon(n), rot),
        "star": lambda n, seed: rotate_points(star_polygon(n), rot),
        "random": lambda n, seed: rotate_points(random_polygon(n, seed), rot),
        # NOTE: "spiral" is intentionally omitted from the default suite because many
        # triangulation baselines interpret this generator as non-simple on some sizes.
    }

    # Keep sizes modest by default; can be overridden via --sizes.
    sizes = list(args.sizes)
    polygons_per_config = int(args.polygons_per_config)

    timeouts = {
        "ours": float(args.timeout_ours),
        "garey": float(args.timeout_garey),
        "hertel": float(args.timeout_hertel),
        "seidel": float(args.timeout_seidel),
    }

    rows = []
    summary = {}

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        poly_path = td / "poly.poly"
        out_path = td / "out.tri"

        total_instances = len(polygon_types) * len(sizes) * polygons_per_config
        done_instances = 0

        for ptype, gen in polygon_types.items():
            for n in sizes:
                times: Dict[str, List[float]] = {k: [] for k in required.keys()}
                reflex_vals: List[int] = []

                for seed in range(polygons_per_config):
                    done_instances += 1
                    if done_instances == 1 or done_instances % 10 == 0 or done_instances == total_instances:
                        print(
                            f"[{done_instances:>5}/{total_instances}] {ptype} n={n} seed={seed}",
                            flush=True,
                        )
                    pts = gen(n, seed)
                    write_poly(pts, poly_path)

                    # Run ours first to get reflex count on the instance.
                    ours = run_cli(required["ours"], poly_path, out_path, timeout_s=timeouts["ours"])
                    if ours.triangles != n - 2:
                        raise RuntimeError(f"ours produced {ours.triangles} triangles for n={n}")
                    if ours.reflex_count is not None:
                        reflex_vals.append(ours.reflex_count)
                    times["ours"].append(ours.time_ms)

                    for alg in ["garey", "hertel", "seidel"]:
                        if alg not in required:
                            continue
                        rr = safe_run_cli(required[alg], poly_path, out_path, timeout_s=timeouts[alg])
                        if rr is None:
                            continue
                        if rr.triangles != n - 2:
                            # Treat as failed run; do not poison the whole benchmark.
                            continue
                        times[alg].append(rr.time_ms)

                r_avg = statistics.mean(reflex_vals) if reflex_vals else 0.0
                record = {"type": ptype, "n": n, "r_avg": r_avg, "r_ratio": (r_avg / n if n else 0.0)}
                for alg, vals in times.items():
                    if vals:
                        record[f"{alg}_mean"] = statistics.mean(vals)
                        record[f"{alg}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
                    else:
                        record[f"{alg}_mean"] = None
                        record[f"{alg}_std"] = None
                summary[(ptype, n)] = record

    # Write results.txt
    lines: List[str] = []
    lines.append("POLYGON TRIANGULATION BENCHMARK (paper)\n")
    lines.append(f"polygons_per_config={polygons_per_config}\n")
    lines.append("Algorithms:\n")
    lines.append("- ours: reflex_cli\n")
    lines.append("- garey: polypartition_mono_cli (Triangulate_MONO)\n")
    lines.append("- hertel: polypartition_hm_cli (ConvexPartition_HM + triangulate pieces)\n")
    if not args.skip_seidel:
        lines.append("- seidel: seidel_cli (palmerc/Seidel)\n")
    lines.append("\n")

    def pm3(mean: Optional[float], std: Optional[float]) -> str:
        if mean is None:
            return "--"
        if std is None:
            return f"{mean:.3f}"
        return f"{mean:.3f}±{std:.3f}"

    alg_order = ["ours", "garey", "hertel"]
    if not args.skip_seidel:
        alg_order.append("seidel")

    for ptype in ["convex", "dent", "random", "star"]:
        lines.append(f"{ptype.upper()}:\n")
        lines.append("n\t r(avg)\t " + "\t ".join(f"{a}(ms)" for a in alg_order) + "\n")
        for n in sizes:
            rec = summary[(ptype, n)]
            cells = [pm3(rec.get(f"{a}_mean"), rec.get(f"{a}_std")) for a in alg_order]
            lines.append(f"{n}\t {int(rec['r_avg'])}\t " + "\t ".join(cells) + "\n")
        lines.append("\n")

    OUT_TXT.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_TXT}")
    write_tex_tables(summary, sizes, ["convex", "dent", "random", "star"])
    print(f"Wrote {OUT_DIR / 'benchmark_table.tex'}")
    print(f"Wrote {OUT_DIR / 'benchmark_full.tex'}")


if __name__ == "__main__":
    main()

