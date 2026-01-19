#!/bin/bash
cd "$(dirname "$0")"

# Generate a simple polygon
python3 << 'EOF'
import random, math

def simple_random_poly(n, seed=0):
    rng = random.Random(seed + n * 1000)
    angles = sorted([rng.random() * 2 * math.pi for _ in range(n)])
    pts = []
    for angle in angles:
        r = 50 + rng.random() * 50
        pts.append((r * math.cos(angle), r * math.sin(angle)))
    return pts

pts = simple_random_poly(50, 0)
with open("/tmp/test.poly", "w") as f:
    f.write(f"{len(pts)}\n")
    for x, y in pts:
        f.write(f"{x} {y}\n")
print("Generated /tmp/test.poly")
EOF

echo "Running chain_optimal_cli..."
./code/ours/chain_optimal_cli --input /tmp/test.poly --output /tmp/test.tri
echo "Exit code: $?"

echo ""
echo "Output file:"
head -20 /tmp/test.tri
