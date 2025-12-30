from note.triangulation import *

points = create_example_polygon()
print("Original points (should be CW):")
for i, p in enumerate(points):
    print(f"  v{i}: {p}")

area = signed_area(points)
print(f"\nSigned area: {area:.2f} ({'CCW' if area > 0 else 'CW'})")

tri = PolygonTriangulator(points)
dec = tri.decomposer
dec.classify_vertices()

print("\nAfter classification (points are now CCW due to reversal):")
for i, v in enumerate(dec.V):
    print(f"  idx {i}: ({v.x}, {v.y}) - {v.vtype.name}")

print(f"\nReflex vertices: {dec.reflex_indices}")

# Run decomposition
dec.decompose_into_chains()
dec.compute_left_edges()
dec.build_helper_sequences()
dec.compute_predecessors()

print(f"\nLeft edges:")
for v_idx in sorted(dec.left_edge.keys()):
    print(f"  Vertex {v_idx}: left edge = {dec.left_edge[v_idx]}")

print(f"\nHelper sequences H(e):")
for edge, seq in dec.helper_seq.items():
    types = [dec.V[i].vtype.name for i in seq]
    print(f"  Edge {edge}: {seq}")
    print(f"    Types: {types}")

print(f"\nPredecessors:")
for v_idx in sorted(dec.predecessor.keys()):
    pred = dec.predecessor[v_idx]
    v = dec.V[v_idx]
    pv = dec.V[pred]
    print(f"  Vertex {v_idx} ({v.vtype.name}): pred = {pred} ({pv.vtype.name})")

dec.compute_diagonal_targets()
print(f"\nDiagonals: {dec.diagonals}")

print("\nExpected from paper figure:")
print("  Original v9 (split) -> idx 2 in reversed -> connects to predecessor")
print("  Original v4 (merge) -> idx 7 in reversed -> connects if pred is merge")  
print("  Original v2 (merge) -> idx 9 in reversed -> connects if pred is merge")
print("  Paper shows diagonals: v9-v2 and v2-v4")
print("  In reversed indices: (2, 9) and (9, 7)")
