#!/usr/bin/env python3
"""
Accurate benchmark using CLI-reported times (no I/O overhead).
Generates multiple random polygons per configuration.
"""
import subprocess
import math
import random
import statistics
import tempfile
from pathlib import Path

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
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_star(n, seed):
    """Generate star polygon with alternating radii."""
    pts = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        r = 100 if i % 2 == 0 else 40
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def gen_spiral(n, seed):
    """Generate spiral polygon (monotone, r=0)."""
    pts = []
    for i in range(n):
        t = i / n
        angle = 4 * math.pi * t
        r = 20 + 80 * t
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts


def write_polygon(pts, path):
    with open(path, 'w') as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")


def run_algo(cli_name, input_path, output_path):
    """Run algorithm and extract time from CLI output."""
    cli = BIN_DIR / cli_name
    if not cli.exists():
        return None, None
    
    try:
        result = subprocess.run(
            [str(cli), '-i', str(input_path), '-o', str(output_path)],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            output = result.stdout
            time_ms = None
            r = 0
            
            for line in output.split('\n'):
                if 'time_ms=' in line:
                    time_ms = float(line.split('time_ms=')[1].split(',')[0].split()[0])
                if 'reflex_count=' in line:
                    r = int(line.split('reflex_count=')[1].split(',')[0].split()[0])
            
            return time_ms, r
    except:
        pass
    return None, None


def main():
    print("=" * 120)
    print("ACCURATE BENCHMARK - Using CLI-reported times (excludes I/O overhead)")
    print("=" * 120)
    
    polygon_types = {
        'convex': gen_convex,
        'spiral': gen_spiral,
        'random': gen_random_star,
        'star': gen_star,
    }
    
    sizes = [1000, 5000, 10000, 20000, 50000, 100000]
    polygons_per_config = 30  # 30 random polygons per configuration
    
    results = {}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        input_path = tmpdir / "input.poly"
        output_path = tmpdir / "output.tri"
        
        for ptype, gen_func in polygon_types.items():
            print(f"\n{ptype.upper()} POLYGONS ({polygons_per_config} polygons per size):")
            print("-" * 110)
            
            for n in sizes:
                ours_times = []
                earcut_times = []
                r_values = []
                
                for seed in range(polygons_per_config):
                    pts = gen_func(n, seed)
                    write_polygon(pts, input_path)
                    
                    # Our algorithm
                    t_ours, r = run_algo('reflex_cli', input_path, output_path)
                    if t_ours is not None:
                        ours_times.append(t_ours)
                        r_values.append(r)
                    
                    # Earcut
                    t_earcut, _ = run_algo('earcut_cli', input_path, output_path)
                    if t_earcut is not None:
                        earcut_times.append(t_earcut)
                
                if ours_times and earcut_times:
                    avg_r = statistics.mean(r_values)
                    avg_ours = statistics.mean(ours_times)
                    std_ours = statistics.stdev(ours_times) if len(ours_times) > 1 else 0
                    avg_earcut = statistics.mean(earcut_times)
                    std_earcut = statistics.stdev(earcut_times) if len(earcut_times) > 1 else 0
                    
                    speedup = avg_earcut / avg_ours if avg_ours > 0 else 0
                    
                    results[(ptype, n)] = {
                        'r_avg': avg_r,
                        'r_ratio': avg_r / n,
                        'ours_mean': avg_ours,
                        'ours_std': std_ours,
                        'earcut_mean': avg_earcut,
                        'earcut_std': std_earcut,
                        'speedup': speedup,
                        'samples': len(ours_times),
                    }
                    
                    print(f"  n={n:>6}  r_avg={int(avg_r):>6} ({100*avg_r/n:>5.1f}%)  "
                          f"Ours={avg_ours:>8.2f}±{std_ours:>6.2f}ms  "
                          f"Earcut={avg_earcut:>8.2f}±{std_earcut:>6.2f}ms  "
                          f"Speedup={speedup:>6.2f}x")
    
    # Generate LaTeX tables
    print("\n" + "=" * 120)
    print("LATEX TABLE FOR MAIN SECTION")
    print("=" * 120)
    
    print("""
\\begin{table}[t]
\\centering
\\caption{Running times (ms) averaged over 30 random polygons per configuration.}
\\label{tab:results}
\\smallskip
\\begin{tabular}{lrrrrrrr}
\\toprule
\\textbf{Type} & $n$ & $r$ (avg) & $r/n$ & \\textbf{Ours} & \\textbf{Earcut} & \\textbf{Speedup} \\\\
\\midrule""")
    
    for ptype in ['convex', 'spiral', 'random', 'star']:
        for n in [10000, 50000, 100000]:
            key = (ptype, n)
            if key in results:
                r = results[key]
                print(f"{ptype.title()} & {n:,} & {int(r['r_avg']):,} & {100*r['r_ratio']:.0f}\\% & "
                      f"${r['ours_mean']:.1f} \\pm {r['ours_std']:.1f}$ & "
                      f"${r['earcut_mean']:.1f} \\pm {r['earcut_std']:.1f}$ & "
                      f"{r['speedup']:.1f}$\\times$ \\\\")
    
    print("""\\bottomrule
\\end{tabular}
\\end{table}""")
    
    # Appendix table
    print("\n" + "=" * 120)
    print("LATEX TABLE FOR APPENDIX (All sizes)")
    print("=" * 120)
    
    print("""
\\begin{table}[h]
\\centering
\\caption{Extended benchmark: running times (ms) averaged over 30 polygons. Standard deviations in parentheses.}
\\label{tab:extended}
\\smallskip
\\small
\\begin{tabular}{llrrrrrr}
\\toprule
\\textbf{Type} & $n$ & $r$ (avg) & $r/n$ & \\textbf{Ours (ms)} & \\textbf{Earcut (ms)} & \\textbf{Speedup} \\\\
\\midrule""")
    
    for ptype in ['convex', 'spiral', 'random', 'star']:
        for n in sizes:
            key = (ptype, n)
            if key in results:
                r = results[key]
                print(f"{ptype.title()} & {n:,} & {int(r['r_avg']):,} & {100*r['r_ratio']:.0f}\\% & "
                      f"{r['ours_mean']:.2f} ({r['ours_std']:.2f}) & "
                      f"{r['earcut_mean']:.2f} ({r['earcut_std']:.2f}) & "
                      f"{r['speedup']:.1f}$\\times$ \\\\")
    
    print("""\\bottomrule
\\end{tabular}
\\end{table}""")


if __name__ == '__main__':
    main()

