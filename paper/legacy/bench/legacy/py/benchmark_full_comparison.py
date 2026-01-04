#!/usr/bin/env python3
"""
Full comparison benchmark with multiple random polygons per configuration.
Compares: Ours, Garey, Hertel, Kirkpatrick, Earcut
"""
import subprocess
import time
import math
import random
import statistics
from pathlib import Path

from garey_impl import triangulate_garey, triangulate_hertel, triangulate_kirkpatrick

ROOT = Path(__file__).parent.parent
BIN_DIR = ROOT / 'build' / 'bin'


def gen_random_star(n, seed):
    """Generate random star-shaped polygon."""
    random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = random.uniform(30, 100)
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_convex(n, seed):
    """Generate convex polygon."""
    random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_star(n, seed):
    """Generate star polygon with alternating radii."""
    random.seed(seed)
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100 if i % 2 == 0 else 40
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_spiral(n, seed):
    """Generate spiral polygon (monotone, r=0)."""
    random.seed(seed)
    pts = []
    for i in range(n):
        t = i / n
        angle = 4 * math.pi * t
        r = 20 + 80 * t
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def write_polygon(pts, path):
    """Write polygon to file."""
    with open(path, 'w') as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")


def run_our_cpp(pts, tmpdir):
    """Run our C++ implementation."""
    input_path = tmpdir / "input.poly"
    output_path = tmpdir / "output.tri"
    write_polygon(pts, input_path)
    
    try:
        t0 = time.perf_counter()
        result = subprocess.run(
            [str(BIN_DIR / 'reflex_cli'), '-i', str(input_path), '-o', str(output_path)],
            capture_output=True, text=True, timeout=60
        )
        elapsed = (time.perf_counter() - t0) * 1000
        
        if result.returncode == 0:
            output = result.stdout
            r = 0
            if 'reflex_count=' in output:
                r = int(output.split('reflex_count=')[1].split(',')[0])
            return elapsed, r
    except:
        pass
    return None, None


def run_earcut_cpp(pts, tmpdir):
    """Run Earcut C++ implementation."""
    input_path = tmpdir / "input.poly"
    output_path = tmpdir / "output.tri"
    write_polygon(pts, input_path)
    
    try:
        t0 = time.perf_counter()
        result = subprocess.run(
            [str(BIN_DIR / 'earcut_cli'), '-i', str(input_path), '-o', str(output_path)],
            capture_output=True, text=True, timeout=60
        )
        elapsed = (time.perf_counter() - t0) * 1000
        
        if result.returncode == 0:
            return elapsed
    except:
        pass
    return None


def run_python_algo(algo_func, pts):
    """Run Python baseline algorithm."""
    t0 = time.perf_counter()
    triangles, r = algo_func(pts)
    elapsed = (time.perf_counter() - t0) * 1000
    return elapsed, r, len(triangles)


def main():
    import tempfile
    
    print("=" * 110)
    print("COMPREHENSIVE BENCHMARK: Ours vs Garey vs Hertel vs Kirkpatrick vs Earcut")
    print("=" * 110)
    
    # Configuration
    polygon_types = {
        'convex': gen_convex,
        'spiral': gen_spiral,
        'random': gen_random_star,
        'star': gen_star,
    }
    
    sizes = [1000, 5000, 10000, 20000, 50000]
    polygons_per_config = 20  # Generate 20 different random polygons per config
    
    results = {}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        for ptype, gen_func in polygon_types.items():
            print(f"\n{ptype.upper()} POLYGONS ({polygons_per_config} polygons per size):")
            print("-" * 100)
            
            for n in sizes:
                # Collect times across multiple random polygons
                ours_times = []
                earcut_times = []
                garey_times = []
                hertel_times = []
                kirk_times = []
                r_values = []
                
                for seed in range(polygons_per_config):
                    pts = gen_func(n, seed)
                    
                    # Our algorithm (C++)
                    t_ours, r = run_our_cpp(pts, tmpdir)
                    if t_ours is not None:
                        ours_times.append(t_ours)
                        r_values.append(r)
                    
                    # Earcut (C++)
                    t_earcut = run_earcut_cpp(pts, tmpdir)
                    if t_earcut is not None:
                        earcut_times.append(t_earcut)
                    
                    # Python baselines (only for smaller sizes due to speed)
                    if n <= 10000:
                        t_garey, _, _ = run_python_algo(triangulate_garey, list(pts))
                        garey_times.append(t_garey)
                        
                        t_hertel, _, _ = run_python_algo(triangulate_hertel, list(pts))
                        hertel_times.append(t_hertel)
                        
                        t_kirk, _, _ = run_python_algo(triangulate_kirkpatrick, list(pts))
                        kirk_times.append(t_kirk)
                
                # Compute averages
                avg_r = statistics.mean(r_values) if r_values else 0
                avg_ours = statistics.mean(ours_times) if ours_times else float('nan')
                avg_earcut = statistics.mean(earcut_times) if earcut_times else float('nan')
                avg_garey = statistics.mean(garey_times) if garey_times else float('nan')
                avg_hertel = statistics.mean(hertel_times) if hertel_times else float('nan')
                avg_kirk = statistics.mean(kirk_times) if kirk_times else float('nan')
                
                std_ours = statistics.stdev(ours_times) if len(ours_times) > 1 else 0
                std_earcut = statistics.stdev(earcut_times) if len(earcut_times) > 1 else 0
                
                # Store results
                results[(ptype, n)] = {
                    'r_avg': avg_r,
                    'r_ratio': avg_r / n,
                    'ours': avg_ours,
                    'ours_std': std_ours,
                    'earcut': avg_earcut,
                    'earcut_std': std_earcut,
                    'garey': avg_garey,
                    'hertel': avg_hertel,
                    'kirk': avg_kirk,
                    'samples': len(ours_times),
                }
                
                # Compute speedups
                speedup_earcut = avg_earcut / avg_ours if avg_ours > 0 else 0
                speedup_garey = avg_garey / avg_ours if avg_ours > 0 and not math.isnan(avg_garey) else 0
                
                print(f"  n={n:>6}  r_avg={int(avg_r):>5} ({100*avg_r/n:>4.1f}%)  "
                      f"Ours={avg_ours:>8.2f}ms  Earcut={avg_earcut:>8.2f}ms ({speedup_earcut:>5.1f}x)  "
                      f"Garey={avg_garey:>8.2f}ms ({speedup_garey:>5.1f}x)")
    
    # Generate LaTeX tables
    print("\n" + "=" * 110)
    print("LATEX TABLE FOR MAIN SECTION")
    print("=" * 110)
    
    print("""
\\begin{table}[t]
\\centering
\\caption{Running times (ms) averaged over 20 random polygons per configuration.}
\\label{tab:results}
\\smallskip
\\begin{tabular}{lrrrrrr}
\\toprule
\\textbf{Type} & $n$ & $r$ (avg) & $r/n$ & \\textbf{Ours} & \\textbf{Earcut} & \\textbf{Speedup} \\\\
\\midrule""")
    
    for ptype in ['convex', 'spiral', 'random', 'star']:
        for n in [10000, 50000]:
            key = (ptype, n)
            if key in results:
                r = results[key]
                speedup = r['earcut'] / r['ours'] if r['ours'] > 0 else 0
                print(f"{ptype.title()} & {n:,} & {int(r['r_avg']):,} & {100*r['r_ratio']:.0f}\\% & "
                      f"{r['ours']:.1f} & {r['earcut']:.1f} & {speedup:.1f}$\\times$ \\\\")
    
    print("""\\bottomrule
\\end{tabular}
\\end{table}""")
    
    # Generate appendix table with more detail
    print("\n" + "=" * 110)
    print("LATEX TABLE FOR APPENDIX (Extended)")
    print("=" * 110)
    
    print("""
\\begin{table}[h]
\\centering
\\caption{Extended benchmark: running times (ms) averaged over 20 polygons per configuration.}
\\label{tab:extended}
\\smallskip
\\small
\\begin{tabular}{llrrrrrrr}
\\toprule
\\textbf{Type} & $n$ & $r$ & $r/n$ & \\textbf{Ours} & \\textbf{Earcut} & \\textbf{Garey} & \\textbf{Hertel} & \\textbf{Kirk.} \\\\
\\midrule""")
    
    for ptype in ['convex', 'spiral', 'random', 'star']:
        for n in sizes:
            key = (ptype, n)
            if key in results:
                r = results[key]
                garey_str = f"{r['garey']:.1f}" if not math.isnan(r['garey']) else "---"
                hertel_str = f"{r['hertel']:.1f}" if not math.isnan(r['hertel']) else "---"
                kirk_str = f"{r['kirk']:.1f}" if not math.isnan(r['kirk']) else "---"
                
                print(f"{ptype.title()} & {n:,} & {int(r['r_avg']):,} & {100*r['r_ratio']:.0f}\\% & "
                      f"{r['ours']:.1f} & {r['earcut']:.1f} & {garey_str} & {hertel_str} & {kirk_str} \\\\")
    
    print("""\\bottomrule
\\end{tabular}
\\end{table}""")
    
    # Summary
    print("\n" + "=" * 110)
    print("SUMMARY")
    print("=" * 110)
    
    for ptype in ['convex', 'spiral', 'random', 'star']:
        for n in [10000, 50000]:
            key = (ptype, n)
            if key in results:
                r = results[key]
                speedup = r['earcut'] / r['ours'] if r['ours'] > 0 else 0
                print(f"{ptype:12} n={n:>6}: r_avg={int(r['r_avg']):>5} ({100*r['r_ratio']:>4.1f}%)  "
                      f"Ours={r['ours']:>6.1f}ms  Speedup over Earcut: {speedup:.1f}x")


if __name__ == '__main__':
    main()

