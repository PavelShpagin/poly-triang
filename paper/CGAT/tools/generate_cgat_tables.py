#!/usr/bin/env python3
"""
Generate CGAT LaTeX tables from the raw benchmark CSV produced by
`scripts/benchmark_paper.py --out-raw-csv ...`.

Expected input CSV schema:
  polygon_type,n,seed,k_count,algorithm,time_ms,triangles,reflex_count

Output:
  - benchmark_table.tex (random)
  - benchmark_full.tex (convex + dent + random + star)
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ALGO_ORDER = ["ours", "seidel", "garey", "hertel"]
ALGO_LABEL = {
    "ours": r"\textbf{Ours}",
    "seidel": "Seidel",
    "garey": "Garey",
    "hertel": "Hertel--Mehlhorn",
}

FAMILY_ORDER = ["convex", "dent", "random", "star"]
FAMILY_LABEL = {
    "convex": "Convex",
    "dent": "Dent",
    "random": "Random",
    "star": "Star",
}

SIZES = [500, 1000, 2000, 5000, 10000]


@dataclass(frozen=True)
class MeanStd:
    mean: float
    std: float


def mean_std(values: List[float]) -> MeanStd:
    if not values:
        raise ValueError("mean_std called with empty list")
    m = statistics.mean(values)
    s = statistics.stdev(values) if len(values) > 1 else 0.0
    return MeanStd(mean=m, std=s)


def fmt_k(ms: MeanStd) -> str:
    return f"{int(round(ms.mean))} $\\pm$ {int(round(ms.std))}"


def fmt_time(ms: MeanStd) -> str:
    return f"{ms.mean:.2f} $\\pm$ {ms.std:.2f}"


def bold(s: str, is_min: bool) -> str:
    return f"\\textbf{{{s}}}" if is_min else s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True, help="Raw CSV from benchmark_paper.py")
    ap.add_argument("--outdir", type=Path, required=True, help="Directory to write generated .tex tables")
    args = ap.parse_args()

    # k_count per (ptype, n, seed)
    k_by_seed: Dict[Tuple[str, int, int], int] = {}

    # time_ms per (ptype, n, alg)
    times: Dict[Tuple[str, int, str], List[float]] = {}

    with open(args.input, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            ptype = (r.get("polygon_type") or "").strip()
            if not ptype:
                continue
            n = int(r["n"])
            seed = int(r["seed"])
            alg = (r.get("algorithm") or "").strip()

            k_count = int(r["k_count"]) if (r.get("k_count") or "").strip() else 0
            kkey = (ptype, n, seed)
            prev = k_by_seed.get(kkey)
            if prev is None:
                k_by_seed[kkey] = k_count
            elif prev != k_count:
                raise RuntimeError(f"Inconsistent k_count for {kkey}: {prev} vs {k_count}")

            t_str = (r.get("time_ms") or "").strip()
            if t_str:
                times.setdefault((ptype, n, alg), []).append(float(t_str))

    # Aggregate k by (ptype, n)
    k_by_conf: Dict[Tuple[str, int], MeanStd] = {}
    for (ptype, n, _seed), k in k_by_seed.items():
        k_by_conf.setdefault((ptype, n), []).append(float(k))
    for key, vals in list(k_by_conf.items()):
        k_by_conf[key] = mean_std(vals)  # type: ignore[assignment]

    # Aggregate times by (ptype, n, alg)
    t_by_conf_alg: Dict[Tuple[str, int, str], MeanStd] = {}
    for key, vals in times.items():
        if vals:
            t_by_conf_alg[key] = mean_std(vals)

    args.outdir.mkdir(parents=True, exist_ok=True)

    # -------------------------
    # Main table (Convex + Random) -- formatted like the full table
    # -------------------------
    lines: List[str] = []
    lines.append("\\begin{table*}[t]")
    lines.append("\\centering")
    lines.append("\\caption{Running time (ms) on convex and random polygons (mean $\\pm$ stdev over instances).}")
    lines.append("\\label{tab:benchmark}")
    lines.append("\\small")
    lines.append("\\begin{tabular}{llrrrrr}")
    lines.append("\\toprule")
    lines.append("Type & $n$ & $k$ & \\textbf{Ours} & Seidel & Garey & Hertel--Mehlhorn \\\\")
    lines.append("\\midrule")

    main_families = ["convex", "random"]
    for fi, fam in enumerate(main_families):
        fam_label = FAMILY_LABEL.get(fam, fam.capitalize())
        first = True
        for n in SIZES:
            kstats = k_by_conf.get((fam, n))
            k_cell = fmt_k(kstats) if kstats else "--"

            means: Dict[str, float] = {}
            cells: Dict[str, str] = {}
            for alg in ALGO_ORDER:
                st = t_by_conf_alg.get((fam, n, alg))
                if st is None:
                    cells[alg] = "--"
                    continue
                means[alg] = st.mean
                cells[alg] = fmt_time(st)

            mn = min(means.values()) if means else math.inf
            for alg in ALGO_ORDER:
                if alg in means:
                    cells[alg] = bold(cells[alg], abs(means[alg] - mn) <= 1e-12)

            type_cell = fam_label if first else ""
            lines.append(
                f"{type_cell} & {n:,} & {k_cell} & {cells['ours']} & {cells['seidel']} & {cells['garey']} & {cells['hertel']} \\\\"
            )
            first = False
        if fi + 1 < len(main_families):
            lines.append("\\midrule")

    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("\\end{table*}")
    (args.outdir / "benchmark_table.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # -------------------------
    # Full table (appendix)
    # -------------------------
    full: List[str] = []
    # Use a non-floating placement in the appendix to avoid large whitespace
    # between the appendix heading and the table.
    full.append("\\begin{table}[H]")
    full.append("\\centering")
    full.append("\\caption{Running time (ms) across polygon families (mean over instances).}")
    full.append("\\label{tab:benchmark-full}")
    full.append("\\small")
    full.append("\\begin{tabular}{llrrrrr}")
    full.append("\\toprule")
    full.append("Type & $n$ & $k$ & \\textbf{Ours} & Seidel & Garey & Hertel--Mehlhorn \\\\")
    full.append("\\midrule")

    for fam_i, fam in enumerate(FAMILY_ORDER):
        first = True
        for n in SIZES:
            kstats = k_by_conf.get((fam, n))
            k_cell = fmt_k(kstats) if kstats else "--"

            means: Dict[str, float] = {}
            cells: Dict[str, str] = {}
            for alg in ALGO_ORDER:
                st = t_by_conf_alg.get((fam, n, alg))
                if st is None:
                    cells[alg] = "--"
                    continue
                means[alg] = st.mean
                cells[alg] = fmt_time(st)

            mn = min(means.values()) if means else math.inf
            for alg in ALGO_ORDER:
                if alg in means:
                    cells[alg] = bold(cells[alg], abs(means[alg] - mn) <= 1e-12)

            type_cell = FAMILY_LABEL.get(fam, fam.capitalize()) if first else ""
            full.append(
                f"{type_cell} & {n:,} & {k_cell} & {cells['ours']} & {cells['seidel']} & {cells['garey']} & {cells['hertel']} \\\\"
            )
            first = False
        if fam_i + 1 < len(FAMILY_ORDER):
            full.append("\\midrule")

    full.append("\\bottomrule")
    full.append("\\end{tabular}")
    full.append("\\end{table}")
    (args.outdir / "benchmark_full.tex").write_text("\n".join(full) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

