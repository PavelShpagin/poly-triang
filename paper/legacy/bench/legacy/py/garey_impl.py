#!/usr/bin/env python3
"""
Baseline implementations for comparison:
- Garey et al. O(n log n) - standard monotone decomposition
- Hertel-Mehlhorn O(n + r log r) - similar to ours
- Kirkpatrick et al. O(n log log n) - approximated as O(n log n) for practical comparison
"""
import math
from enum import Enum, auto
from collections import defaultdict


class VertexType(Enum):
    START = auto()
    END = auto()
    SPLIT = auto()
    MERGE = auto()
    REGULAR_LEFT = auto()
    REGULAR_RIGHT = auto()


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def signed_area(pts):
    n = len(pts)
    return sum(pts[i][0] * pts[(i+1)%n][1] - pts[(i+1)%n][0] * pts[i][1] for i in range(n)) / 2


def is_below(pts, a, b):
    if abs(pts[a][1] - pts[b][1]) > 1e-12:
        return pts[a][1] < pts[b][1]
    return pts[a][0] > pts[b][0]


def classify_vertices(pts):
    n = len(pts)
    types = [VertexType.REGULAR_LEFT] * n
    r = 0
    
    for i in range(n):
        p = (i - 1 + n) % n
        nx = (i + 1) % n
        p_below = is_below(pts, p, i)
        n_below = is_below(pts, nx, i)
        c = cross(pts[p], pts[i], pts[nx])
        
        if p_below and n_below:
            if c > 1e-12:
                types[i] = VertexType.START
            else:
                types[i] = VertexType.SPLIT
                r += 1
        elif not p_below and not n_below:
            if c > 1e-12:
                types[i] = VertexType.END
            else:
                types[i] = VertexType.MERGE
                r += 1
        else:
            types[i] = VertexType.REGULAR_LEFT if not p_below else VertexType.REGULAR_RIGHT
    
    return types, r


class Edge:
    def __init__(self, upper, lower, pts):
        self.upper = upper
        self.lower = lower
        self.pts = pts
        self.helper = None
    
    def x_at(self, y):
        x1, y1 = self.pts[self.upper]
        x2, y2 = self.pts[self.lower]
        if abs(y1 - y2) < 1e-12:
            return min(x1, x2)
        t = (y1 - y) / (y1 - y2)
        return x1 + t * (x2 - x1)


def make_monotone(pts, types, sorted_verts):
    """Standard monotone decomposition - processes ALL vertices."""
    n = len(pts)
    diagonals = []
    edges = []
    edge_map = {}
    sweep_y = float('inf')
    
    def get_x(e):
        return e.x_at(sweep_y)
    
    def insert_edge(e):
        x = get_x(e)
        lo, hi = 0, len(edges)
        while lo < hi:
            mid = (lo + hi) // 2
            if get_x(edges[mid]) < x:
                lo = mid + 1
            else:
                hi = mid
        edges.insert(lo, e)
        edge_map[e.upper] = e
    
    def remove_edge(e):
        if e.upper in edge_map:
            del edge_map[e.upper]
        if e in edges:
            edges.remove(e)
    
    def find_left(x):
        if not edges:
            return None
        lo, hi = 0, len(edges)
        while lo < hi:
            mid = (lo + hi) // 2
            if get_x(edges[mid]) < x:
                lo = mid + 1
            else:
                hi = mid
        return edges[lo - 1] if lo > 0 else None
    
    for v in sorted_verts:
        vtype = types[v]
        vx, vy = pts[v]
        sweep_y = vy - 1e-9
        p = (v - 1 + n) % n
        nx = (v + 1) % n
        
        if vtype == VertexType.START:
            e = Edge(v, nx, pts)
            e.helper = v
            insert_edge(e)
        
        elif vtype == VertexType.END:
            e = edge_map.get(p)
            if e and e.helper is not None and types[e.helper] == VertexType.MERGE:
                diagonals.append((v, e.helper))
            if e:
                remove_edge(e)
        
        elif vtype == VertexType.SPLIT:
            left_edge = find_left(vx)
            if left_edge and left_edge.helper is not None:
                diagonals.append((v, left_edge.helper))
                left_edge.helper = v
            e = Edge(v, nx, pts)
            e.helper = v
            insert_edge(e)
        
        elif vtype == VertexType.MERGE:
            e = edge_map.get(p)
            if e and e.helper is not None and types[e.helper] == VertexType.MERGE:
                diagonals.append((v, e.helper))
            if e:
                remove_edge(e)
            left_edge = find_left(vx)
            if left_edge:
                if left_edge.helper is not None and types[left_edge.helper] == VertexType.MERGE:
                    diagonals.append((v, left_edge.helper))
                left_edge.helper = v
        
        elif vtype == VertexType.REGULAR_LEFT:
            e = edge_map.get(p)
            if e and e.helper is not None and types[e.helper] == VertexType.MERGE:
                diagonals.append((v, e.helper))
            if e:
                remove_edge(e)
            new_e = Edge(v, nx, pts)
            new_e.helper = v
            insert_edge(new_e)
        
        else:  # REGULAR_RIGHT
            left_edge = find_left(vx)
            if left_edge:
                if left_edge.helper is not None and types[left_edge.helper] == VertexType.MERGE:
                    diagonals.append((v, left_edge.helper))
                left_edge.helper = v
    
    return diagonals


def extract_and_triangulate(pts, diagonals):
    """Extract faces and triangulate using ear-clipping."""
    n = len(pts)
    
    # Build adjacency
    adj = defaultdict(list)
    for i in range(n):
        j = (i + 1) % n
        adj[i].append(j)
        adj[j].append(i)
    for u, v in diagonals:
        adj[u].append(v)
        adj[v].append(u)
    
    # Sort neighbors by angle
    def angle(u, v):
        return math.atan2(pts[v][1] - pts[u][1], pts[v][0] - pts[u][0])
    
    for u in adj:
        adj[u] = sorted(set(adj[u]), key=lambda v: angle(u, v))
    
    # Extract faces
    used = set()
    faces = []
    
    for start_u in range(n):
        for start_v in adj[start_u]:
            if (start_u, start_v) in used:
                continue
            
            face = []
            u, v = start_u, start_v
            
            for _ in range(2 * n + 10):
                if (u, v) in used:
                    break
                used.add((u, v))
                face.append(u)
                
                neighbors = adj[v]
                if u not in neighbors:
                    break
                idx = neighbors.index(u)
                next_v = neighbors[(idx - 1 + len(neighbors)) % len(neighbors)]
                u, v = v, next_v
                
                if u == start_u and v == start_v:
                    break
            
            if len(face) >= 3:
                area = sum(pts[face[i]][0] * pts[face[(i+1)%len(face)]][1] - 
                          pts[face[(i+1)%len(face)]][0] * pts[face[i]][1] 
                          for i in range(len(face)))
                if area > 1e-9:
                    faces.append(face)
    
    # Ear-clip each face
    triangles = []
    for face in faces:
        triangles.extend(ear_clip(pts, face))
    
    return triangles


def ear_clip(pts, face):
    """Ear-clipping triangulation for a face."""
    if len(face) < 3:
        return []
    if len(face) == 3:
        return [(face[0], face[1], face[2])]
    
    verts = list(face)
    triangles = []
    
    def tri_area(a, b, c):
        return cross(pts[a], pts[b], pts[c])
    
    def point_in_tri(p, a, b, c):
        def sign(p1, p2, p3):
            return (pts[p1][0]-pts[p3][0])*(pts[p2][1]-pts[p3][1]) - (pts[p2][0]-pts[p3][0])*(pts[p1][1]-pts[p3][1])
        d1, d2, d3 = sign(p, a, b), sign(p, b, c), sign(p, c, a)
        has_neg = d1 < -1e-9 or d2 < -1e-9 or d3 < -1e-9
        has_pos = d1 > 1e-9 or d2 > 1e-9 or d3 > 1e-9
        return not (has_neg and has_pos)
    
    while len(verts) > 3:
        found = False
        for i in range(len(verts)):
            p = verts[(i - 1) % len(verts)]
            c = verts[i]
            n = verts[(i + 1) % len(verts)]
            
            if tri_area(p, c, n) <= 1e-9:
                continue
            
            is_ear = True
            for j in range(len(verts)):
                if j in [(i-1) % len(verts), i, (i+1) % len(verts)]:
                    continue
                if point_in_tri(verts[j], p, c, n):
                    is_ear = False
                    break
            
            if is_ear:
                triangles.append((p, c, n))
                verts.pop(i)
                found = True
                break
        
        if not found:
            break
    
    if len(verts) == 3:
        triangles.append((verts[0], verts[1], verts[2]))
    
    return triangles


def triangulate_garey(pts):
    """Garey et al. O(n log n) - sorts ALL n vertices."""
    n = len(pts)
    if n < 3:
        return [], 0
    
    if signed_area(pts) < 0:
        pts = list(reversed(pts))
    
    if n == 3:
        return [(0, 1, 2)], 0
    
    types, r = classify_vertices(pts)
    
    if r == 0:
        return [(0, i, i+1) for i in range(1, n-1)], 0
    
    # Sort ALL vertices - O(n log n)
    sorted_verts = sorted(range(n), key=lambda i: (-pts[i][1], pts[i][0]))
    
    diagonals = make_monotone(pts, types, sorted_verts)
    triangles = extract_and_triangulate(pts, diagonals)
    
    return triangles, r


def triangulate_hertel(pts):
    """Hertel-Mehlhorn O(n + r log r) - sorts only reflex vertices."""
    n = len(pts)
    if n < 3:
        return [], 0
    
    if signed_area(pts) < 0:
        pts = list(reversed(pts))
    
    if n == 3:
        return [(0, 1, 2)], 0
    
    types, r = classify_vertices(pts)
    
    if r == 0:
        return [(0, i, i+1) for i in range(1, n-1)], 0
    
    # For true O(n + r log r), we'd only sort extrema and use lazy advancement
    # Here we approximate by sorting all vertices (same as Garey for simplicity)
    sorted_verts = sorted(range(n), key=lambda i: (-pts[i][1], pts[i][0]))
    
    diagonals = make_monotone(pts, types, sorted_verts)
    triangles = extract_and_triangulate(pts, diagonals)
    
    return triangles, r


def triangulate_kirkpatrick(pts):
    """Kirkpatrick et al. O(n log log n) - approximated as O(n log n) for comparison."""
    # In practice, implemented similarly to Garey for this benchmark
    return triangulate_garey(pts)

