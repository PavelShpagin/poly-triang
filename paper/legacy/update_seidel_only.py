#!/usr/bin/env python3
"""
Update paper tables by rerunning ONLY the Seidel baseline on the paper benchmark suite.

This script:
- Reads existing `paper/results.txt` for (ours, garey, hertel) means/stds.
- Regenerates the *same* polygon instances (convex/dent/random/star; sizes from results.txt;
  seeds = [0..polygons_per_config-1]) using the generators in `paper/benchmark.py`.
- Runs `build/bin/seidel_cli` and computes mean/std per (type, n).
- Rewrites:
  - `paper/results.txt` (drops Earcut, adds Seidel column)
  - `paper/generated/benchmark_table.tex` (random-only table; drops Earcut)
  - `paper/generated/benchmark_full.tex` (all families table; drops Earcut)
"""

from __future__ import annotations

import re
import statistics
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


ROOT = Path(__file__).resolve().parent.parent
PAPER_DIR = Path(__file__).resolve().parent
BIN_DIR = ROOT / "build" / "bin"

IN_RESULTS = PAPER_DIR / "results.txt"
OUT_RESULTS = PAPER_DIR / "results.txt"
OUT_DIR = PAPER_DIR / "generated"

SEIDEL_EXE = BIN_DIR / "seidel_cli"


@dataclass
class PM:
    mean: Optional[float]
    std: Optional[float]


def _parse_pm(s: str) -> PM:
    s = s.strip()
    if s == "--":
        return PM(None, None)
    # supports: "1.234±0.056" or "1.234 ± 0.056"
    if "±" in s:
        a, b = s.split("±", 1)
        return PM(float(a.strip()), float(b.strip()))
    return PM(float(s), 0.0)


def _fmt_pm_txt(pm: PM) -> str:
    if pm.mean is None:
        return "--"
    if pm.std is None:
        return f"{pm.mean:.3f}"
    return f"{pm.mean:.3f}±{pm.std:.3f}"


def _fmt_pm_tex(pm: PM) -> str:
    if pm.mean is None:
        return "--"
    if pm.std is None:
        return f"{pm.mean:.2f}"
    return f"{pm.mean:.2f} $\\pm$ {pm.std:.2f}"


def _read_results_without_recomputing(path: Path) -> Tuple[int, List[int], Dict[Tuple[str, int], dict]]:
    """
    Parse the existing results.txt that may include an Earcut column, and return:
    - polygons_per_config
    - sizes (sorted as they appear in the file)
    - summary[(ptype, n)] with keys:
        r_avg, ours_mean/std, garey_mean/std, hertel_mean/std
    """
    txt = path.read_text(encoding="utf-8")
    m = re.search(r"polygons_per_config=(\d+)", txt)
    if not m:
        raise RuntimeError(f"Could not parse polygons_per_config from {path}")
    ppc = int(m.group(1))

    summary: Dict[Tuple[str, int], dict] = {}
    sizes_seen: List[int] = []

    lines = [ln.rstrip("\n") for ln in txt.splitlines()]
    i = 0
    current: Optional[str] = None
    while i < len(lines):
        ln = lines[i].strip()
        if ln in ("CONVEX:", "DENT:", "RANDOM:", "STAR:"):
            current = ln[:-1].lower()
            i += 1
            # header line
            if i >= len(lines):
                break
            header = lines[i].strip()
            cols = [c.strip() for c in header.split("\t") if c.strip()]
            # Expect: n, r(avg), ours, [earcut], garey, hertel
            # We'll locate by name prefix.
            idx = {name: j for j, name in enumerate(cols)}
            if "n" not in idx or "r(avg)" not in idx:
                raise RuntimeError(f"Unexpected header in {current}: {header}")
            i += 1
            while i < len(lines):
                row = lines[i].strip()
                if not row:
                    break
                parts = [p.strip() for p in row.split("\t")]
                n = int(parts[idx["n"]])
                r_avg = int(parts[idx["r(avg)"]])
                ours = _parse_pm(parts[idx.get("ours(ms)", idx.get("ours", -1))])
                garey = _parse_pm(parts[idx.get("garey(ms)", idx.get("garey", -1))])
                hertel = _parse_pm(parts[idx.get("hertel(ms)", idx.get("hertel", -1))])
                summary[(current, n)] = {
                    "type": current,
                    "n": n,
                    "r_avg": r_avg,
                    "ours": ours,
                    "garey": garey,
                    "hertel": hertel,
                }
                if n not in sizes_seen:
                    sizes_seen.append(n)
                i += 1
            current = None
        i += 1

    if not summary:
        raise RuntimeError(f"No benchmark rows parsed from {path}")

    return ppc, sizes_seen, summary


def _write_tables(summary: Dict[Tuple[str, int], dict], sizes: List[int]) -> None:
    OUT_DIR.mkdir(exist_ok=True)

    # Random-only table.
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
        cols = [("ours", rec["ours"]), ("seidel", rec["seidel"]), ("garey", rec["garey"]), ("hertel", rec["hertel"])]
        best = min((pm.mean for _, pm in cols if pm.mean is not None), default=None)
        cells = []
        for _, pm in cols:
            cell = _fmt_pm_tex(pm)
            if best is not None and pm.mean is not None and abs(pm.mean - best) < 1e-12:
                cell = r"\textbf{" + cell + "}"
            cells.append(cell)
        lines.append(f"{n:,} & {int(rec['r_avg'])} & " + " & ".join(cells) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    (OUT_DIR / "benchmark_table.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Full table.
    lines = []
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{Running time (ms) across polygon families (mean $\pm$ stdev over instances).}")
    lines.append(r"\label{tab:benchmark-full}")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{llrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"Type & $n$ & $r$ & \textbf{Ours} & Seidel & Garey & Hertel--Mehlhorn \\")
    lines.append(r"\midrule")

    for ptype in ["convex", "dent", "random", "star"]:
        first = True
        for n in sizes:
            rec = summary.get((ptype, n))
            if not rec:
                continue
            cols = [("ours", rec["ours"]), ("seidel", rec["seidel"]), ("garey", rec["garey"]), ("hertel", rec["hertel"])]
            best = min((pm.mean for _, pm in cols if pm.mean is not None), default=None)
            vals = []
            for _, pm in cols:
                cell = _fmt_pm_tex(pm)
                if best is not None and pm.mean is not None and abs(pm.mean - best) < 1e-12:
                    cell = r"\textbf{" + cell + "}"
                vals.append(cell)
            type_col = ptype.capitalize() if first else ""
            first = False
            lines.append(f"{type_col} & {n:,} & {int(rec['r_avg'])} & " + " & ".join(vals) + r" \\")
        lines.append(r"\midrule")
    if lines[-1] == r"\midrule":
        lines[-1] = r"\bottomrule"
    else:
        lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table*}")
    (OUT_DIR / "benchmark_full.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_results_txt(ppc: int, sizes: List[int], summary: Dict[Tuple[str, int], dict]) -> None:
    lines: List[str] = []
    lines.append("POLYGON TRIANGULATION BENCHMARK (paper)\n")
    lines.append(f"polygons_per_config={ppc}\n")
    lines.append("Algorithms:\n")
    lines.append("- ours: reflex_cli\n")
    lines.append("- garey: polypartition_mono_cli (Triangulate_MONO)\n")
    lines.append("- hertel: polypartition_hm_cli (ConvexPartition_HM + triangulate pieces)\n")
    lines.append("- seidel: seidel_cli (palmerc/Seidel)\n")
    lines.append("\n")

    for ptype in ["convex", "dent", "random", "star"]:
        lines.append(f"{ptype.upper()}:\n")
        lines.append("n\t r(avg)\t ours(ms)\t seidel(ms)\t garey(ms)\t hertel(ms)\n")
        for n in sizes:
            rec = summary.get((ptype, n))
            if not rec:
                continue
            lines.append(
                f"{n}\t {int(rec['r_avg'])}\t "
                + "\t ".join(
                    [
                        _fmt_pm_txt(rec["ours"]),
                        _fmt_pm_txt(rec["seidel"]),
                        _fmt_pm_txt(rec["garey"]),
                        _fmt_pm_txt(rec["hertel"]),
                    ]
                )
                + "\n"
            )
        lines.append("\n")

    OUT_RESULTS.write_text("".join(lines), encoding="utf-8")


def main() -> None:
    if not IN_RESULTS.exists():
        raise SystemExit(f"Missing {IN_RESULTS}. Run the paper benchmark first.")
    if not SEIDEL_EXE.exists():
        raise SystemExit(f"Missing {SEIDEL_EXE}. Build first.")

    # Parse baseline numbers (ours/garey/hertel) from current results.txt.
    ppc, sizes, summary = _read_results_without_recomputing(IN_RESULTS)

    # Import polygon generators/helpers from paper/benchmark.py.
    import benchmark  # type: ignore

    rot = 0.123456789
    polygon_types = {
        "convex": lambda n, seed: benchmark.rotate_points(benchmark.convex_polygon(n), rot),
        "dent": lambda n, seed: benchmark.rotate_points(benchmark.dent_polygon(n), rot),
        "star": lambda n, seed: benchmark.rotate_points(benchmark.star_polygon(n), rot),
        "random": lambda n, seed: benchmark.rotate_points(benchmark.random_polygon(n, seed), rot),
    }

    # Rerun ONLY seidel.
    print(f"Running Seidel only: {len(polygon_types)} types × {len(sizes)} sizes × {ppc} seeds = {len(polygon_types)*len(sizes)*ppc} runs")
    # Keep this tight: we only want Seidel numbers when it's responsive.
    timeout_s = 1.0

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        poly_path = td / "poly.poly"
        out_path = td / "out.tri"

        total = len(polygon_types) * len(sizes) * ppc
        done = 0

        for ptype, gen in polygon_types.items():
            for n in sizes:
                vals: List[float] = []
                for seed in range(ppc):
                    done += 1
                    if done == 1 or done % 10 == 0 or done == total:
                        print(f"[{done:>5}/{total}] seidel {ptype} n={n} seed={seed}", flush=True)
                    pts = gen(n, seed)
                    benchmark.write_poly(pts, poly_path)
                    rr = benchmark.safe_run_cli(SEIDEL_EXE, poly_path, out_path, timeout_s=timeout_s)
                    if rr is None:
                        continue
                    if rr.triangles != n - 2:
                        continue
                    vals.append(rr.time_ms)

                if vals:
                    mean = statistics.mean(vals)
                    std = statistics.stdev(vals) if len(vals) > 1 else 0.0
                    pm = PM(mean, std)
                else:
                    pm = PM(None, None)

                key = (ptype, n)
                if key not in summary:
                    summary[key] = {"type": ptype, "n": n, "r_avg": 0, "ours": PM(None, None), "garey": PM(None, None), "hertel": PM(None, None)}
                summary[key]["seidel"] = pm

    # Ensure all records have a seidel field (may be missing if results.txt didn't include a row).
    for key, rec in summary.items():
        rec.setdefault("seidel", PM(None, None))

    # Rewrite artifacts (drop Earcut; add Seidel).
    _write_results_txt(ppc, sizes, summary)
    _write_tables(summary, sizes)
    print(f"Wrote {OUT_RESULTS}")
    print(f"Wrote {OUT_DIR / 'benchmark_table.tex'}")
    print(f"Wrote {OUT_DIR / 'benchmark_full.tex'}")


if __name__ == "__main__":
    main()


