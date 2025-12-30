
import math
import random
import time
from typing import List, Tuple, Dict, Set

class Vertex:
    def __init__(self, x, y, i):
        self.x = x
        self.y = y
        self.i = i
        self.prev = None
        self.next = None
        self.is_reflex = False

def cross(o, a, b):
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)

class GridEarClipping:
    def __init__(self, points):
        self.n = len(points)
        self.verts = [Vertex(x, y, i) for i, (x, y) in enumerate(points)]
        self.reflex_verts = []
        self.ears = []
        
        # Build Linked List
        for i in range(self.n):
            v = self.verts[i]
            v.prev = self.verts[(i-1)%self.n]
            v.next = self.verts[(i+1)%self.n]
            
        # Classify and Grid
        self.min_x = min(p[0] for p in points)
        self.max_x = max(p[0] for p in points)
        self.min_y = min(p[1] for p in points)
        self.max_y = max(p[1] for p in points)
        
        # Grid parameters
        # Heuristic: grid size approx sqrt(N)
        grid_dim = int(math.sqrt(self.n))
        if grid_dim < 2: grid_dim = 2
        self.cell_w = (self.max_x - self.min_x) / grid_dim
        self.cell_h = (self.max_y - self.min_y) / grid_dim
        if self.cell_w == 0: self.cell_w = 1.0
        if self.cell_h == 0: self.cell_h = 1.0
        
        self.grid = {} # (ix, iy) -> [Vertex]
        
        for v in self.verts:
            if cross(v.prev, v, v.next) < 0:
                v.is_reflex = True
                self.reflex_verts.append(v)
                self.add_to_grid(v)
        
        self.stats = {'checks': 0, 'ears_found': 0}

    def get_cell(self, v):
        ix = int((v.x - self.min_x) / self.cell_w)
        iy = int((v.y - self.min_y) / self.cell_h)
        return (ix, iy)

    def add_to_grid(self, v):
        k = self.get_cell(v)
        if k not in self.grid: self.grid[k] = []
        self.grid[k].append(v)

    def remove_from_grid(self, v):
        k = self.get_cell(v)
        if k in self.grid:
            try:
                self.grid[k].remove(v)
            except ValueError:
                pass 

    def is_ear(self, v):
        if v.is_reflex: return False
        
        # Triangle box
        prev, nxt = v.prev, v.next
        x_min = min(v.x, prev.x, nxt.x)
        x_max = max(v.x, prev.x, nxt.x)
        y_min = min(v.y, prev.y, nxt.y)
        y_max = max(v.y, prev.y, nxt.y)
        
        # Check grid cells intersecting box
        ix_min = int((x_min - self.min_x) / self.cell_w)
        ix_max = int((x_max - self.min_x) / self.cell_w)
        iy_min = int((y_min - self.min_y) / self.cell_h)
        iy_max = int((y_max - self.min_y) / self.cell_h)
        
        # Bounding box of triangle
        
        for ix in range(ix_min, ix_max + 1):
            for iy in range(iy_min, iy_max + 1):
                if (ix, iy) in self.grid:
                    for r in self.grid[(ix, iy)]:
                        if r == v or r == prev or r == nxt: continue
                        self.stats['checks'] += 1
                        # Check if r in triangle
                        # Barycentric or simply 3 cross products
                        # Point in Triangle Test
                        c1 = cross(prev, v, r)
                        c2 = cross(v, nxt, r)
                        c3 = cross(nxt, prev, r)
                        # Assuming CCW triangle (convex vertex)
                        # All must be >= 0 (left or on)
                        if c1 >= 0 and c2 >= 0 and c3 >= 0:
                            return False
        return True

    def triangulate(self):
        diagonals = []
        # Initial Ear Scan
        # List of candidate ears?
        # A simple queue
        queue = [v for v in self.verts if not v.is_reflex]
        
        # To avoid O(N^2), we process queue.
        # But if we process v and it's not an ear, we might re-check it later?
        # Standard Ear Clipping:
        # Maintain a doubly linked list.
        # Check current v. If ear, clip, check neighbors.
        # If not ear, move to next.
        
        # For O(N) behavior, we need to visit each vertex O(1) times.
        # If we skip v, when do we revisit?
        # Only when a neighbor is clipped.
        
        curr = self.verts[0]
        count = self.n
        
        # FIST strategy:
        # Maintain "CHECK_NEXT" flag?
        
        # Simple Loop
        # While n > 3:
        #   check curr
        #   if ear: clip, update neighbors, stay at neighbor?
        #   else: next
        
        # Iteration limit to prevent infinite loop on bad geometry
        iters = 0
        limit = self.n * self.n 
        
        while count > 3 and iters < limit:
            iters += 1
            if self.is_ear(curr):
                # Clip
                prev, nxt = curr.prev, curr.next
                diagonals.append((prev.i, nxt.i))
                
                # Remove curr
                prev.next = nxt
                nxt.prev = prev
                
                # Update reflex status of neighbors
                # If neighbor was reflex, it might become convex
                if prev.is_reflex:
                    if cross(prev.prev, prev, prev.next) >= 0:
                        prev.is_reflex = False
                        self.remove_from_grid(prev)
                if nxt.is_reflex:
                    if cross(nxt.prev, nxt, nxt.next) >= 0:
                        nxt.is_reflex = False
                        self.remove_from_grid(nxt)
                        
                curr = nxt # Continue from neighbor
                count -= 1
            else:
                curr = curr.next
                
        return diagonals, self.stats

# Benchmark Code
def create_comb_jitter(n, seed=42):
    random.seed(seed)
    teeth = n // 4
    pts = []
    pts.append((-10, -50))
    for i in range(teeth):
        x_base = i * 4
        pts.append((x_base, random.random() * 20)) 
        pts.append((x_base + 1, 150 + random.random())) 
        pts.append((x_base + 2, random.random() * 20)) 
        pts.append((x_base + 3, random.random() * 10)) 
    pts.append((teeth*4 + 10, -50))
    return pts

def benchmark():
    print("=" * 80)
    print("GRID EAR CLIPPING BENCHMARK")
    print("=" * 80)
    print(f"{'n':>6} {'Checks':>10} {'Chk/N':>8} {'Time(s)':>10}")
    
    for n in [1000, 5000, 10000, 20000]:
        pts = create_comb_jitter(n)
        t0 = time.time()
        ear = GridEarClipping(pts)
        _, stats = ear.triangulate()
        t1 = time.time()
        
        checks = stats['checks']
        print(f"{n:>6} {checks:>10} {checks/n:>8.2f} {t1-t0:>10.4f}")

if __name__ == "__main__":
    benchmark()

