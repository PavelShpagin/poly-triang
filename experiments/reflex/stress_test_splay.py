
import math
import random
import time
from triangulate_splay import SplayTriangulator

def create_comb_jitter(n, seed=42):
    """
    A comb with perturbed y-coordinates to force non-sequential access 
    in the sweep line status structure.
    """
    random.seed(seed)
    pts = []
    teeth = n // 4
    
    # We want many vertical edges active.
    # Top vertices at y ~ 100
    # Bottom vertices at y ~ 0, but perturbed
    
    # Generating a comb-like structure
    # Up, Down, Right, Up, Down, Right...
    
    width = 1000.0
    dx = width / teeth
    
    tops = []
    bottoms = []
    
    for i in range(teeth):
        x_left = i * dx
        x_right = x_left + dx * 0.5
        
        # Bottom-Left
        y_bl = random.random() * 10.0
        # Top-Left
        y_tl = 100.0 + random.random() * 10.0
        # Top-Right
        y_tr = 100.0 + random.random() * 10.0
        # Bottom-Right
        y_br = random.random() * 10.0
        
        # We need to order them to form a polygon.
        # Let's just make a simple zig-zag path:
        # (0,0) -> (0,100) -> (1,100) -> (1,0) -> (2,0) -> ...
        
        # But we want the "bottom" turns to be processed in random order of x.
        # So we assign random y's to the turns.
        
    # Better construction:
    # 1. Create upper chain (monotonic x)
    # 2. Create lower chain (monotonic x)
    # 3. Connect them.
    # But this doesn't create "vertical" bars that stay in the sweep line.
    
    # Back to Comb:
    # Vertices: 
    # 0: (0, y_0_low)
    # 1: (0, y_0_high)
    # 2: (1, y_0_high)
    # 3: (1, y_0_low)
    # 4: (2, y_1_low) ...
    
    xs = []
    ys = []
    
    for i in range(teeth):
        x = i * 2
        # Up leg
        pts.append((x, 100 + random.random()))
        pts.append((x+0.5, 100 + random.random()))
        
    # Now go back down
    # This is not a comb, this is a strip.
    # A comb needs to go up and down.
    
    pts = []
    # Left wall
    pts.append((-10, 0))
    pts.append((-10, 200))
    
    # Top wall
    pts.append((teeth * 4 + 10, 200))
    
    # Right wall
    pts.append((teeth * 4 + 10, 0))
    
    # Bottom "Teeth"
    # We want N/2 vertices here.
    # Going from right to left? Or left to right?
    # Polygon order must be CCW.
    # So we are at bottom-right, going left.
    
    bottom_pts = []
    for i in range(teeth):
        x = i * 4
        # We want a spike going up.
        # (x, 0) -> (x+1, 100) -> (x+2, 0)
        # But we want the "tips" and "bases" to be at random heights 
        # so the sweep line hits them in random x-order?
        
        # If we have vertical bars, they are active for a long range of y.
        # The events are the vertices.
        # If we have 1000 bars, we have 1000 edges in the BST.
        # If the vertices at the bottom are at y=0.1, y=0.5, y=0.2...
        # We process them in y-order.
        # y=0.1 is at x=50
        # y=0.2 is at x=900
        # y=0.3 is at x=10
        # This forces the BST search to jump around.
        
        # Bases:
        y_base_l = random.random() * 10
        y_high = 100 + random.random() * 10
        y_base_r = random.random() * 10
        
        # To maintain simple polygon, the "teeth" shouldn't intersect.
        # x is increasing.
        # We are going Right to Left.
        
        # Let's generate Left to Right and reverse.
        
        # Tooth i:
        # p1 = (x, y_base_l)
        # p2 = (x+1, y_high)
        # p3 = (x+2, y_base_r)
        
        bottom_pts.append((x, y_base_l))
        bottom_pts.append((x+1, y_high))
        bottom_pts.append((x+2, y_base_r))
        # spacer
        bottom_pts.append((x+3, random.random()*5))
        
    # Reverse bottom points to go Right to Left
    # Actually, let's just reverse the logic to append correctly.
    # Or just use the points we generated and fix the hull.
    
    # Let's just make a "Comb" where the bottom vertices have random Ys.
    # And the top is a flat line.
    
    final_pts = []
    # Start bottom left
    final_pts.append((-10, -50))
    
    # Teeth
    for i in range(teeth):
        x_base = i * 4
        # Up
        final_pts.append((x_base, random.random() * 20)) # Low
        final_pts.append((x_base + 1, 150 + random.random())) # High
        final_pts.append((x_base + 2, random.random() * 20)) # Low
        # Spacer
        final_pts.append((x_base + 3, random.random() * 10)) # Low
        
    # End bottom right
    final_pts.append((teeth*4 + 10, -50))
    
    # Close the loop
    # This is CCW?
    # x increases. y varies.
    # Then back to start.
    # This might cross if we are not careful.
    # The "Low" points are around y=0-20.
    # The "High" points are around y=150.
    # The return path is at y=-50.
    # Safe.
    
    # Is it CCW?
    # Low->High->Low->Spacer... moving right.
    # Then down to -50.
    # Then left to start.
    # Then up to first Low.
    # Yes, CCW (interior is to the left).
    # Wait, if we walk right along bottom, interior is Up. Left is Up. Correct.
    
    return final_pts

def benchmark_jitter():
    print("=" * 80)
    print("SPLAY KILLER BENCHMARK (JITTER COMB)")
    print("=" * 80)
    print(f"{'n':>6} {'ops':>8} {'ops/n':>8}")
    
    for n in [100, 500, 1000, 2000, 5000, 10000]:
        pts = create_comb_jitter(n)
        tri = SplayTriangulator(pts)
        _, stats = tri.triangulate()
        ops = stats['splay_ops']
        print(f"{len(pts):>6} {ops:>8} {ops/len(pts):>8.2f}")

if __name__ == "__main__":
    benchmark_jitter()

