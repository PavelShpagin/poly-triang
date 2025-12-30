#!/usr/bin/env python3
"""
Statistical benchmark with confidence intervals.
"""
import subprocess
import csv
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).parent.parent
BIN_DIR = ROOT / 'build' / 'bin'
POLY_DIR = ROOT / 'polygons' / 'generated'
RESULTS_DIR = ROOT / 'results'

def run_algo(algo, input_path, output_path):
    """Run algorithm and return time in ms."""
    cli = BIN_DIR / f"{algo}_cli"
    if not cli.exists():
        return None, None
    
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
                    return time_ms, r
    except subprocess.TimeoutExpired:
        pass
    return None, None

def main():
    print("=" * 100)
    print("STATISTICAL BENCHMARK WITH CONFIDENCE INTERVALS")
    print("=" * 100)
    
    RESULTS_DIR.mkdir(exist_ok=True)
    
    # Focus on key polygon types
    polygon_types = ['convex', 'spiral', 'random', 'star']
    sizes = [10000, 50000, 100000]
    algorithms = ['reflex', 'earcut']
    runs_per_config = 20  # For statistical stability
    
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
                    t, r = run_algo(algo, poly_path, output_path)
                    if t is not None:
                        times[algo].append(t)
                        if algo == 'reflex':
                            r_val = r
            
            if len(times['reflex']) >= 5 and len(times['earcut']) >= 5:
                # Calculate statistics
                reflex_mean = statistics.mean(times['reflex'])
                reflex_std = statistics.stdev(times['reflex'])
                earcut_mean = statistics.mean(times['earcut'])
                earcut_std = statistics.stdev(times['earcut'])
                
                speedup = earcut_mean / reflex_mean if reflex_mean > 0 else 0
                
                # 95% confidence interval for speedup (approximate)
                reflex_se = reflex_std / math.sqrt(len(times['reflex']))
                earcut_se = earcut_std / math.sqrt(len(times['earcut']))
                
                results[(ptype, n)] = {
                    'reflex_mean': reflex_mean,
                    'reflex_std': reflex_std,
                    'earcut_mean': earcut_mean,
                    'earcut_std': earcut_std,
                    'r': r_val,
                    'speedup': speedup,
                    'samples': len(times['reflex'])
                }
                
                print(f"  n={n:>7}  r={r_val:>6} ({100*r_val/n:>4.1f}%)")
                print(f"    Ours:   {reflex_mean:>8.2f} +/- {reflex_std:>6.2f} ms")
                print(f"    Earcut: {earcut_mean:>8.2f} +/- {earcut_std:>6.2f} ms")
                print(f"    Speedup: {speedup:>6.2f}x ({len(times['reflex'])} samples)")
    
    # Generate final summary
    print("\n" + "=" * 100)
    print("FINAL SUMMARY")
    print("=" * 100)
    
    print(f"\nKey Results at n = 100,000:")
    for ptype in polygon_types:
        key = (ptype, 100000)
        if key in results:
            r = results[key]
            print(f"  {ptype.title():12}: {r['speedup']:.1f}x speedup (r/n = {100*r['r']/100000:.0f}%)")
    
    # Save to CSV
    csv_path = RESULTS_DIR / 'statistical_benchmark.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['polygon_type', 'n', 'r', 'r_ratio', 'reflex_mean', 'reflex_std', 
                        'earcut_mean', 'earcut_std', 'speedup', 'samples'])
        for (ptype, n), r in results.items():
            writer.writerow([ptype, n, r['r'], r['r']/n, r['reflex_mean'], r['reflex_std'],
                           r['earcut_mean'], r['earcut_std'], r['speedup'], r['samples']])
    
    print(f"\nResults saved to {csv_path}")

if __name__ == '__main__':
    main()

