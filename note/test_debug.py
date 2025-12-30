#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')
from note.triangulation import PolygonTriangulator, create_paper_example, VType

pts = create_paper_example()
pt = PolygonTriangulator(pts)
pt._classify_vertices()

def prv(i): return (i - 1 + pt.n) % pt.n
def nxt(i): return (i + 1) % pt.n

print("=== Vertex Classification ===")
for v in pt.V:
    print(f"v{v.idx}: ({v.x:.1f}, {v.y:.1f}) {v.vtype.name}")

print("\n=== Edge Analysis ===")
# Edge i goes from vertex i to vertex nxt(i)
# It is "active" when sweep line is between the y-coords of its endpoints
for i in range(pt.n):
    j = nxt(i)
    y_i, y_j = pt.pts[i][1], pt.pts[j][1]
    if y_i > y_j:
        print(f"Edge {i}: v{i}->v{j}, upper=v{i}, active from y={y_i:.1f} to y={y_j:.1f}")
    else:
        print(f"Edge {i}: v{i}->v{j}, upper=v{j}, active from y={y_j:.1f} to y={y_i:.1f}")

print("\n=== Monotone Decomposition Trace ===")
order = sorted(range(pt.n), key=lambda i: (-pt.pts[i][1], pt.pts[i][0]))

active_edges = []
helper = {}
diagonals = []

def edge_x(upper, y):
    lower = nxt(upper)
    x1, y1 = pt.pts[upper]
    x2, y2 = pt.pts[lower]
    if abs(y1 - y2) < 1e-12:
        return min(x1, x2)
    t = (y - y1) / (y2 - y1)
    return x1 + t * (x2 - x1)

def find_left_edge(v_idx, sweep_y):
    vx = pt.pts[v_idx][0]
    best = -1
    best_x = float('-inf')
    for e in active_edges:
        ex = edge_x(e, sweep_y)
        if ex < vx + 1e-9 and ex > best_x:
            best_x = ex
            best = e
    return best

for v in order:
    vtype = pt.V[v].vtype
    sweep_y = pt.pts[v][1] - 1e-9
    
    print(f"\nProcessing v{v} ({pt.pts[v][0]:.1f}, {pt.pts[v][1]:.1f}) {vtype.name}")
    print(f"  Active edges: {active_edges}")
    print(f"  Helpers: {helper}")
    
    if vtype == VType.START:
        active_edges.append(v)
        helper[v] = v
        print(f"  -> Added edge {v}, helper[{v}] = {v}")
        
    elif vtype == VType.END:
        e = prv(v)
        print(f"  -> Edge ending: {e} (from v{e} to v{v})")
        if e in helper and pt.V[helper[e]].vtype == VType.MERGE:
            diagonals.append((v, helper[e]))
            print(f"  -> DIAGONAL: v{v} - v{helper[e]}")
        if e in active_edges:
            active_edges.remove(e)
            print(f"  -> Removed edge {e}")
        else:
            print(f"  -> Edge {e} not in active_edges!")
            
    elif vtype == VType.SPLIT:
        ej = find_left_edge(v, sweep_y)
        print(f"  -> Left edge: {ej}")
        if ej >= 0 and ej in helper:
            diagonals.append((v, helper[ej]))
            print(f"  -> DIAGONAL: v{v} - v{helper[ej]}")
            helper[ej] = v
        active_edges.append(v)
        helper[v] = v
        print(f"  -> Added edge {v}")
        
    elif vtype == VType.MERGE:
        e = prv(v)
        print(f"  -> Edge ending: {e}")
        if e in helper and pt.V[helper[e]].vtype == VType.MERGE:
            diagonals.append((v, helper[e]))
            print(f"  -> DIAGONAL: v{v} - v{helper[e]}")
        if e in active_edges:
            active_edges.remove(e)
        
        ej = find_left_edge(v, sweep_y)
        print(f"  -> Left edge: {ej}")
        if ej >= 0:
            if ej in helper and pt.V[helper[ej]].vtype == VType.MERGE:
                diagonals.append((v, helper[ej]))
                print(f"  -> DIAGONAL: v{v} - v{helper[ej]}")
            helper[ej] = v
            print(f"  -> Updated helper[{ej}] = {v}")
            
    elif vtype == VType.REG_R:
        e = prv(v)
        print(f"  -> Edge ending: {e}")
        if e in helper and pt.V[helper[e]].vtype == VType.MERGE:
            diagonals.append((v, helper[e]))
            print(f"  -> DIAGONAL: v{v} - v{helper[e]}")
        if e in active_edges:
            active_edges.remove(e)
        active_edges.append(v)
        helper[v] = v
        print(f"  -> Replaced edge {e} with edge {v}")
        
    elif vtype == VType.REG_L:
        ej = find_left_edge(v, sweep_y)
        print(f"  -> Left edge: {ej}")
        if ej >= 0:
            if ej in helper and pt.V[helper[ej]].vtype == VType.MERGE:
                diagonals.append((v, helper[ej]))
                print(f"  -> DIAGONAL: v{v} - v{helper[ej]}")
            helper[ej] = v
            print(f"  -> Updated helper[{ej}] = {v}")

print(f"\n=== Final Diagonals: {diagonals} ===")
print(f"Expected: 3 diagonals for 3 reflex vertices")

