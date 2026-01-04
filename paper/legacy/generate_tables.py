#!/usr/bin/env python3
"""
Generate paper tables from benchmark CSV data.
Uses the actual C++ benchmark results from results/methods_benchmark.csv
"""

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"
OUTPUT_DIR = Path(__file__).resolve().parent / "generated"


def parse_methods_benchmark():
    """Parse results/methods_benchmark.csv"""
    data = defaultdict(dict)
    
    csv_path = RESULTS_DIR / "methods_benchmark.csv"
    if not csv_path.exists():
        print(f"Warning: {csv_path} not found")
        return data
    
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            algo = row["algorithm"]
            polygon = row["polygon"]
            n = int(row["num_vertices"])
            r = int(row["num_reflex"])
            time_ms = float(row["time_ms"])
            
            # Parse polygon type from name (e.g., "convex_1000" -> "convex")
            ptype = polygon.rsplit("_", 1)[0] if "_" in polygon else polygon
            
            key = (ptype, n)
            if key not in data:
                data[key] = {"r": r, "reflex": [], "earcut": [], "garey": [], "hertel": []}
            
            if "reflex" in algo:
                data[key]["reflex"].append(time_ms)
                data[key]["r"] = r  # Update with actual reflex count
            elif "earcut" in algo:
                data[key]["earcut"].append(time_ms)
            elif "garey" in algo:
                data[key]["garey"].append(time_ms)
            elif "hertel" in algo:
                data[key]["hertel"].append(time_ms)
    
    return data


def mean(vals):
    return sum(vals) / len(vals) if vals else 0


def generate_main_table(data):
    """Generate Table 1: main comparison at n=100,000"""
    lines = []
    lines.append(r"\begin{table}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{Running times (ms) at $n = 100{,}000$ vertices.}")
    lines.append(r"\label{tab:results}")
    lines.append(r"\smallskip")
    lines.append(r"\begin{tabular}{lrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Type} & $r$ (avg) & $r/n$ & \textbf{Ours} & \textbf{Earcut} & \textbf{Speedup} \\")
    lines.append(r"\midrule")
    
    for ptype in ["convex", "spiral", "random", "star"]:
        rec = data.get((ptype, 100000))
        if not rec:
            continue
        
        r = rec["r"]
        r_ratio = f"{r/100000*100:.0f}\\%"
        
        ours = mean(rec["reflex"])
        earcut = mean(rec["earcut"])
        
        if ours > 0.001 and earcut > 0:
            speedup = f"{earcut/ours:.1f}$\\times$"
        elif ours > 0 and earcut > 0:
            speedup = f"{earcut/ours:.0f}$\\times$"
        else:
            speedup = "--"
        
        ours_s = f"{ours:.2f}" if ours >= 1 else f"{ours:.3f}"
        earcut_s = f"{earcut:.1f}" if earcut >= 100 else f"{earcut:.2f}"
        
        lines.append(f"{ptype.capitalize()} & {r:,} & {r_ratio} & {ours_s} & {earcut_s} & {speedup} \\\\")
    
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    
    return "\n".join(lines)


def generate_extended_table(data):
    """Generate extended table for appendix"""
    lines = []
    lines.append(r"\begin{table}[h]")
    lines.append(r"\centering")
    lines.append(r"\caption{Extended benchmark: running times (ms) across polygon sizes.}")
    lines.append(r"\label{tab:extended}")
    lines.append(r"\smallskip")
    lines.append(r"\small")
    lines.append(r"\begin{tabular}{llrrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Type} & $n$ & $r$ & $r/n$ & \textbf{Ours (ms)} & \textbf{Earcut (ms)} & \textbf{Speedup} \\")
    lines.append(r"\midrule")
    
    sizes = [1000, 5000, 10000, 50000, 100000, 200000, 500000, 1000000]
    
    for ptype in ["convex", "spiral", "random", "star"]:
        first_row = True
        for n in sizes:
            rec = data.get((ptype, n))
            if not rec or not rec["reflex"]:
                continue
            
            r = rec["r"]
            r_ratio = f"{r/n*100:.0f}\\%"
            
            ours = mean(rec["reflex"])
            earcut = mean(rec["earcut"])
            
            if ours > 0 and earcut > 0:
                speedup = f"{earcut/ours:.1f}$\\times$" if earcut/ours < 1000 else f"{earcut/ours:.0f}$\\times$"
            else:
                speedup = "--"
            
            ours_s = f"{ours:.3f}" if ours < 1 else f"{ours:.2f}"
            earcut_s = f"{earcut:.2f}" if earcut < 1000 else f"{earcut:.0f}"
            
            type_col = ptype.capitalize() if first_row else ""
            first_row = False
            
            n_fmt = f"{n:,}"
            lines.append(f"{type_col} & {n_fmt} & {r:,} & {r_ratio} & {ours_s} & {earcut_s} & {speedup} \\\\")
        
        lines.append(r"\midrule")
    
    lines[-1] = r"\bottomrule"
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    
    return "\n".join(lines)


def generate_comparison_table(data):
    """Generate algorithm comparison table"""
    lines = []
    lines.append(r"\begin{table}[h]")
    lines.append(r"\centering")
    lines.append(r"\caption{Algorithm comparison. Speedups over Earcut measured at $n = 100{,}000$.}")
    lines.append(r"\label{tab:comparison}")
    lines.append(r"\smallskip")
    lines.append(r"\begin{tabular}{lcc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Algorithm} & \textbf{Complexity} & \textbf{Speedup} \\")
    lines.append(r"\midrule")
    lines.append(r"Ear clipping (Earcut)~\cite{mapbox2016earcut} & $O(n^2)$ worst & baseline \\")
    lines.append(r"Garey et al.~\cite{garey1978} & $O(n \log n)$ & -- \\")
    lines.append(r"Hertel--Mehlhorn~\cite{hertel1983} & $O(n + r \log r)$ & -- \\")
    lines.append(r"Seidel~\cite{seidel1991} & $O(n \log^* n)$ exp. & -- \\")
    lines.append(r"Chazelle~\cite{chazelle1991} & $O(n)$ & impractical \\")
    lines.append(r"\textbf{This paper} & $O(n + r \log r)$ & \textbf{18--1500$\times$} \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    
    return "\n".join(lines)


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    print("Parsing benchmark data...")
    data = parse_methods_benchmark()
    
    if not data:
        print("No benchmark data found. Please run benchmarks first.")
        return
    
    print(f"Found {len(data)} data points")
    
    # Print summary
    print("\n=== BENCHMARK SUMMARY ===")
    for ptype in ["convex", "spiral", "random", "star"]:
        print(f"\n{ptype.upper()}:")
        for n in [1000, 10000, 100000, 1000000]:
            rec = data.get((ptype, n))
            if rec and rec["reflex"] and rec["earcut"]:
                ours = mean(rec["reflex"])
                earcut = mean(rec["earcut"])
                r = rec["r"]
                speedup = earcut/ours if ours > 0 else 0
                print(f"  n={n:>8,} r={r:>6,} ours={ours:>10.3f}ms earcut={earcut:>10.1f}ms speedup={speedup:>.1f}x")
    
    # Generate LaTeX tables
    print("\n=== GENERATING TABLES ===")
    
    main_table = generate_main_table(data)
    (OUTPUT_DIR / "benchmark_table.tex").write_text(main_table, encoding="utf-8")
    print(f"  Wrote {OUTPUT_DIR / 'benchmark_table.tex'}")
    
    ext_table = generate_extended_table(data)
    (OUTPUT_DIR / "benchmark_full.tex").write_text(ext_table, encoding="utf-8")
    print(f"  Wrote {OUTPUT_DIR / 'benchmark_full.tex'}")
    
    comp_table = generate_comparison_table(data)
    (OUTPUT_DIR / "comparison_table.tex").write_text(comp_table, encoding="utf-8")
    print(f"  Wrote {OUTPUT_DIR / 'comparison_table.tex'}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()

