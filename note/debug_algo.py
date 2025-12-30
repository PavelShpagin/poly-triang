"""Debug the algorithm step by step."""

import sys
sys.path.insert(0, '/home/pavel/dev/poly-triang/note')

from triangulation import PolygonTriangulator, VType

def create_monotone_polygon():
    """Create a simple y-monotone polygon."""
    return [
        (0, 0), (2, 1), (4, 0), (4, 3), (2, 2), (0, 3)
    ]

pts = create_monotone_polygon()
print("Polygon vertices:")
for i, p in enumerate(pts):
    print(f"  {i}: {p}")

pt = PolygonTriangulator(pts)
pt._classify_vertices()

print("\nVertex types:")
for i in range(pt.n):
    print(f"  {i}: {pt.V[i].name}")

print(f"\nExtrema: {pt.extrema}")
print(f"Reflex count r = {pt.r}")

pt._build_chains()
print(f"\nChains ({len(pt.chains)}):")
for cid, chain in pt.chains.items():
    print(f"  Chain {cid}: vertices={chain.vertex_indices}, top_y={chain.top_y}, bot_y={chain.bot_y}")

# Now run the monotone decomposition
pt._make_monotone_optimal()
print(f"\nDiagonals: {pt.diagonals}")

# Extract faces
tris = pt._triangulate_monotone_pieces()
print(f"\nTriangles: {tris}")

