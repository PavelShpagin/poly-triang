#!/usr/bin/env python3
"""
Final benchmark using pre-generated correct polygons.
"""
import subprocess
import sys
import csv
from pathlib import Path

ROOT = Path(__file__).parent.parent
BIN_DIR = ROOT / 'build' / 'bin'
POLY_DIR = ROOT / 'polygons' / 'generated'
RESULTS_DIR = ROOT / 'results'

def run_algo(algo, input_path, output_path):
    """Run algorithm and return time in ms."""
    cli = BIN_DIR / f"{algo}_cli"
    if not cli.exists():
        return None, None, None
    
    try:
        result = subprocess.run(
            [str(cli), '--input', str(input_path), '--output', str(output_path)],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            output = result.stdout
            for line in output.split('\n'):
                if 'time_ms=' in line:
                    time_ms = float(line.split('time_ms=')[1].split(',')[0].split()[0])
                    r = 0
                    if 'reflex_count=' in line:
                        r = int(line.split('reflex_count=')[1].split(',')[0].split()[0])
                    tris = 0
                    if 'triangles=' in line:
                        tris = int(line.split('triangles=')[1].split(',')[0])
                    return time_ms, r, tris
    except subprocess.TimeoutExpired:
        pass
    return None, None, None

def main():
    print("=" * 90)
    print("FINAL BENCHMARK ON PRE-GENERATED POLYGONS")
    print("=" * 90)
    
    RESULTS_DIR.mkdir(exist_ok=True)
    
    polygon_types = ['convex', 'random', 'spiral', 'star']
    sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
    algorithms = ['reflex', 'earcut']
    runs_per_config = 5  # Multiple runs for stability
    
    csv_path = RESULTS_DIR / 'final_correct_benchmark.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['polygon_type', 'n', 'algorithm', 'run', 'time_ms', 'reflex_count', 'triangles'])
    
    results = {}
    
    for ptype in polygon_types:
        print(f"\n{ptype.upper()}:")
        
        for n in sizes:
            poly_path = POLY_DIR / f"{ptype}_{n}.poly"
            if not poly_path.exists():
                continue
            
            output_path = RESULTS_DIR / f"out_{ptype}_{n}.tri"
            
            times = {algo: [] for algo in algorithms}
            r_val = 0
            
            for run in range(runs_per_config):
                for algo in algorithms:
                    t, r, tris = run_algo(algo, poly_path, output_path)
                    if t is not None:
                        times[algo].append(t)
                        if algo == 'reflex':
                            r_val = r
                        
                        with open(csv_path, 'a', newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow([ptype, n, algo, run, t, r or 0, tris or 0])
            
            if times['reflex'] and times['earcut']:
                avg_reflex = sum(times['reflex']) / len(times['reflex'])
                avg_earcut = sum(times['earcut']) / len(times['earcut'])
                speedup = avg_earcut / avg_reflex if avg_reflex > 0 else 0
                
                results[(ptype, n)] = {
                    'reflex': avg_reflex,
                    'earcut': avg_earcut,
                    'r': r_val,
                    'speedup': speedup
                }
                
                print(f"  n={n:>7}  r={r_val:>6}  "
                      f"Ours={avg_reflex:>9.3f}ms  Earcut={avg_earcut:>9.3f}ms  "
                      f"Speedup={speedup:>7.2f}x")
    
    print(f"\nResults saved to {csv_path}")
    
    # Generate LaTeX table
    print("\n" + "=" * 90)
    print("LATEX TABLE FOR MAIN SECTION (Low-reflex polygons)")
    print("=" * 90)
    
    print("\\begin{table}[t]")
    print("\\centering")
    print("\\caption{Running times (ms) on low-reflex polygons.}")
    print("\\label{tab:results}")
    print("\\smallskip")
    print("\\begin{tabular}{lrrrr}")
    print("\\toprule")
    print("\\textbf{Type} & $n$ & \\textbf{Ours} & \\textbf{Earcut} & \\textbf{Speedup} \\\\")
    print("\\midrule")
    
    for ptype in ['convex', 'spiral']:
        for n in [1000, 10000, 50000, 100000]:
            key = (ptype, n)
            if key in results:
                r = results[key]
                print(f"{ptype.title()} & {n:,} & {r['reflex']:.2f} & {r['earcut']:.2f} & {r['speedup']:.1f}$\\times$ \\\\")
    
    print("\\bottomrule")
    print("\\end{tabular}")
    print("\\end{table}")
    
    # Generate appendix table
    print("\n" + "=" * 90)
    print("LATEX TABLE FOR APPENDIX (All polygon types)")
    print("=" * 90)
    
    print("\\begin{table}[h]")
    print("\\centering")
    print("\\caption{Running times (ms) by polygon type.}")
    print("\\label{tab:extended}")
    print("\\smallskip")
    print("\\begin{tabular}{llrrrrr}")
    print("\\toprule")
    print("\\textbf{Type} & $n$ & $r$ & \\textbf{Ours} & \\textbf{Earcut} & \\textbf{Speedup} \\\\")
    print("\\midrule")
    
    for ptype in ['convex', 'spiral', 'random', 'star']:
        for n in [1000, 5000, 10000, 50000]:
            key = (ptype, n)
            if key in results:
                r = results[key]
                print(f"{ptype.title()} & {n:,} & {r['r']:,} & {r['reflex']:.2f} & {r['earcut']:.2f} & {r['speedup']:.1f}$\\times$ \\\\")
    
    print("\\bottomrule")
    print("\\end{tabular}")
    print("\\end{table}")


if __name__ == '__main__':
    main()

