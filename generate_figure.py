import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.patches import FancyBboxPatch

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Polygon vertices (CCW order)
verts = np.array([
    (0, 2.5),    # 0
    (1, 5.5),    # 1
    (2.5, 4),    # 2 - MERGE (both neighbors above, angle > pi)
    (4, 6.5),    # 3
    (5.5, 5),    # 4 - MERGE (both neighbors above, angle > pi)
    (7, 7),      # 5
    (8, 5.5),    # 6
    (6.5, 3.5),  # 7 - regular (one neighbor above, one below)
    (8, 1.5),    # 8
    (5, 2.5),    # 9 - SPLIT (both neighbors below, angle > pi)
    (3, 0),      # 10
    (1.5, 1.5),  # 11 - regular (one neighbor above, one below)
])

n = len(verts)
# Reflex = split + merge vertices only (per algorithm definition)
reflex_idx = {2, 4, 9}
convex_idx = set(range(n)) - reflex_idx

def draw_polygon_base(ax):
    poly = patches.Polygon(verts, facecolor='#f5f5f5', edgecolor='black', linewidth=1.5, zorder=1)
    ax.add_patch(poly)

def draw_vertices(ax, reflex_black=True, all_circles=False):
    for i in range(n):
        x, y = verts[i]
        if all_circles:
            circle = patches.Circle((x, y), 0.18, facecolor='white', edgecolor='black', linewidth=1.2, zorder=10)
            ax.add_patch(circle)
        elif i in reflex_idx:
            size = 0.28
            rect = patches.Rectangle((x - size/2, y - size/2), size, size, 
                                      facecolor='black', edgecolor='black', linewidth=1, zorder=10)
            ax.add_patch(rect)
        else:
            circle = patches.Circle((x, y), 0.18, facecolor='white', edgecolor='black', linewidth=1.2, zorder=10)
            ax.add_patch(circle)

def setup_ax(ax, title):
    ax.set_xlim(-1, 9.5)
    ax.set_ylim(-1, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=12, style='italic', y=-0.08)

# (a) Simple polygon P
ax1 = axes[0]
draw_polygon_base(ax1)
draw_vertices(ax1)
setup_ax(ax1, "(a) Simple polygon $P$")

# (b) Reflex graph G
ax2 = axes[1]
draw_polygon_base(ax2)

# Horizontal chords from reflex vertices (dashed, inside polygon)
# Only at split/merge vertices per algorithm definition
reflex_y = {2: 4, 4: 5, 9: 2.5}

def get_chord_endpoints(y_level):
    intersections = []
    for i in range(n):
        p1 = verts[i]
        p2 = verts[(i + 1) % n]
        if (p1[1] - y_level) * (p2[1] - y_level) < 0:
            t = (y_level - p1[1]) / (p2[1] - p1[1])
            x = p1[0] + t * (p2[0] - p1[0])
            intersections.append(x)
    if len(intersections) >= 2:
        return min(intersections), max(intersections)
    return None, None

for idx, y in reflex_y.items():
    x_left, x_right = get_chord_endpoints(y)
    if x_left is not None:
        ax2.plot([x_left, x_right], [y, y], color='gray', linestyle='--', linewidth=1, zorder=2)

# Reflex graph edges - diagonals output during monotone decomposition
# v2 (merge) finds v4 (merge) as helper of left edge
# v9 (split) outputs diagonal to helper v2
reflex_edges = [
    (2, 4),   # v2 connects to v4 (merge finds merge as helper)
    (9, 2),   # v9 connects to v2 (split outputs diagonal to helper)
]
for i, j in reflex_edges:
    ax2.plot([verts[i][0], verts[j][0]], [verts[i][1], verts[j][1]], 
             color='#555555', linewidth=2, zorder=3)

draw_vertices(ax2)
setup_ax(ax2, "(b) Reflex graph $G$")

# (c) Triangulation of P
ax3 = axes[2]
draw_polygon_base(ax3)

# Proper triangulation using ear clipping / fan decomposition
# First, the monotone decomposition diagonals
monotone_diags = [
    (2, 1),   # from split v3 to helper
    (4, 6),   # from split v5 to helper  
    (7, 4),   # from merge v8
    (9, 7),   # from merge v10
    (11, 9),  # from merge v12
]

# Full triangulation diagonals (9 valid, non-crossing diagonals for 12-vertex polygon)
triangulation = [
    (0, 2),
    (2, 4),
    (5, 7),
    (7, 9),
    (9, 11),
    (2, 11),
    (2, 9),
    (4, 7),
    (2, 7),
]

for i, j in triangulation:
    ax3.plot([verts[i][0], verts[j][0]], [verts[i][1], verts[j][1]], 
             color='#888888', linewidth=1, zorder=2)

draw_vertices(ax3, all_circles=True)
setup_ax(ax3, "(c) Triangulation of $P$")

plt.tight_layout()
plt.savefig('note/figures/algorithm_steps.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
print("Figure saved to note/figures/algorithm_steps.png")
