#!/usr/bin/env python3
"""Generate additional polygon types for benchmark."""
import math
import random
from pathlib import Path

POLY_DIR = Path(__file__).parent.parent / 'polygons' / 'generated'

def gen_comb(n, seed=None):
    """Generate comb-shaped polygon with many teeth (high reflex count)."""
    if seed is not None:
        random.seed(seed)
    
    # Number of teeth
    teeth = n // 4
    pts = []
    
    # Bottom edge left to right
    for i in range(teeth):
        x = i * 4
        pts.append((x, 0))
        pts.append((x + 1, 3))
        pts.append((x + 2, 3))
        pts.append((x + 3, 0))
    
    # Close the polygon
    pts.append((teeth * 4, 0))
    pts.append((teeth * 4, -1))
    pts.append((0, -1))
    
    # Truncate or extend to exactly n vertices
    return pts[:n] if len(pts) >= n else pts + pts[:n - len(pts)]

def gen_zigzag(n, seed=None):
    """Generate zigzag polygon."""
    if seed is not None:
        random.seed(seed)
    
    pts = []
    for i in range(n // 2):
        pts.append((i * 2, 0 if i % 2 == 0 else 2))
    for i in range(n // 2 - 1, -1, -1):
        pts.append((i * 2 + 1, 5 if i % 2 == 0 else 3))
    
    return pts[:n]

def gen_staircase(n, seed=None):
    """Generate staircase polygon."""
    if seed is not None:
        random.seed(seed)
    
    steps = n // 4
    pts = []
    
    # Up stairs
    for i in range(steps):
        pts.append((i * 2, i * 2))
        pts.append((i * 2 + 2, i * 2))
    
    # Top
    pts.append((steps * 2, steps * 2))
    pts.append((steps * 2, steps * 2 + 1))
    
    # Down stairs (shifted)
    for i in range(steps - 1, -1, -1):
        pts.append((i * 2 + 2, i * 2 + 1))
        pts.append((i * 2, i * 2 + 1))
    
    # Bottom
    pts.append((0, 1))
    pts.append((0, 0))
    
    return pts[:n]

def gen_skewed(n, seed=None):
    """Generate skewed polygon (like a parallelogram with jagged edges)."""
    if seed is not None:
        random.seed(seed)
    
    pts = []
    # Bottom edge
    for i in range(n // 2):
        x = i * 2
        y = random.uniform(0, 0.5)
        pts.append((x, y))
    
    # Top edge (reversed, skewed)
    for i in range(n // 2 - 1, -1, -1):
        x = i * 2 + 5  # Skew
        y = 10 + random.uniform(0, 0.5)
        pts.append((x, y))
    
    return pts[:n]

def gen_gear(n, seed=None):
    """Generate gear-shaped polygon with teeth."""
    if seed is not None:
        random.seed(seed)
    
    teeth = n // 4
    pts = []
    
    for i in range(teeth):
        angle = 2 * math.pi * i / teeth
        # Outer point
        pts.append((100 * math.cos(angle), 100 * math.sin(angle)))
        # Inner corner 1
        angle1 = angle + math.pi / (teeth * 2)
        pts.append((60 * math.cos(angle1), 60 * math.sin(angle1)))
        # Inner corner 2
        angle2 = angle + math.pi / teeth - math.pi / (teeth * 2)
        pts.append((60 * math.cos(angle2), 60 * math.sin(angle2)))
        # Next tooth start
        angle3 = angle + math.pi / teeth
        pts.append((100 * math.cos(angle3), 100 * math.sin(angle3)))
    
    return pts[:n]

def write_polygon(pts, path):
    """Write polygon to file."""
    with open(path, 'w') as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x} {y}\n")

def main():
    POLY_DIR.mkdir(exist_ok=True, parents=True)
    
    sizes = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
    generators = {
        'comb': gen_comb,
        'zigzag': gen_zigzag,
        'staircase': gen_staircase,
        'skewed': gen_skewed,
        'gear': gen_gear,
    }
    
    for name, gen in generators.items():
        print(f"Generating {name} polygons...")
        for n in sizes:
            pts = gen(n, seed=42)
            path = POLY_DIR / f"{name}_{n}.poly"
            write_polygon(pts, path)
            print(f"  {path.name}: {len(pts)} vertices")

if __name__ == '__main__':
    main()

