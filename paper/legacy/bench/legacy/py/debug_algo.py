#!/usr/bin/env python3
"""Debug the triangulation algorithm."""
import subprocess
from pathlib import Path

ROOT = Path(__file__).parent.parent
BIN_DIR = ROOT / 'build' / 'bin'
POLY_DIR = ROOT / 'polygons' / 'generated'

# Run our algorithm
for ptype in ['convex', 'spiral', 'random', 'star']:
    poly = POLY_DIR / f"{ptype}_1000.poly"
    result = subprocess.run(
        [str(BIN_DIR / 'reflex_cli'), '-i', str(poly), '-o', '/tmp/test.tri'],
        capture_output=True, text=True
    )
    output = result.stdout.strip()
    
    # Get expected triangles
    with open(poly) as f:
        n = int(f.readline().strip())
    expected = n - 2
    
    # Parse actual
    actual = 0
    if 'triangles=' in output:
        actual = int(output.split('triangles=')[1].split(',')[0])
    
    status = "PASS" if actual == expected else "FAIL"
    print(f"{status}: {ptype} n={n} expected={expected} actual={actual}")

