# Brainstorming: Simpler Linear-Time Triangulation

We seek an $O(n)$ algorithm (or close to it) that is simpler than Chazelle's.
We assume a Real-RAM model or allow practical assumptions (like grids) if justified.

## Idea 1: Splay-Tree Sweep
- **Concept:** Use a Splay Tree as the status structure in a plane sweep.
- **Why:** Splay trees exploit locality. The "working set" of edges in a polygon sweep might be small or structured.
- **Complexity:** Worst case $O(n \log n)$, but potentially $O(n)$ for many polygon classes.
- **Simplicity:** Very high. Just a standard sweep with a splay tree.
- **Theoretical Angle:** Analyze the dynamic finger property in the context of polygon geometry.

## Idea 2: Adaptive Grid / Bucket Sweep
- **Concept:** Overlay a uniform grid on the polygon. Use the grid to index edges.
- **Why:** Ray shooting (finding the edge to the left) becomes $O(1)$ expected if the grid is fine enough.
- **Complexity:** $O(n)$ average. Worst case depends on distribution.
- **Simplicity:** High.
- **Theoretical Angle:** "Fat" polygons or probabilistic analysis.

## Idea 3: Smart Monotone Chain Merging
- **Concept:** Decompose boundary into monotone chains. Merge them intelligently.
- **Why:** Sorting is the bottleneck. If we can avoid full sort or merge efficiently.
- **Complexity:** $O(n \log (\text{chains}))$.
- **Simplicity:** Moderate.

## Idea 4: Reflex-Driven Decomposition (The "Rigorous" Approach)
- **Concept:** Only process reflex vertices. Use horizontal chords.
- **Why:** $r$ is often small.
- **Complexity:** $O(n + r \log r)$.
- **Simplicity:** Moderate. Already explored.

## Idea 5: Topological Sweep of the Dual
- **Concept:** Don't sweep by $y$. Sweep by "visibility depth" or topological ordering.
- **Complexity:** ?
- **Simplicity:** Low.

## Idea 6: Randomized Incremental (Seidel's Variant)
- **Concept:** Insert segments in random order. Build trapezoidal map.
- **Complexity:** $O(n \log^* n)$ expected.
- **Simplicity:** Moderate.
- **Refinement:** Can we derandomize or simplify the data structure?

## Idea 7: Ear-Clipping with Locality
- **Concept:** Maintain a set of ears. Use a "visited" list to avoid re-checking non-ears.
- **Complexity:** Worst case $O(n^2)$ if we pick bad ears.
- **Refinement:** "Safe" ears?

## Idea 8: Vertical Decomposition with "Zone" Theorem
- **Concept:** Insert edges one by one. Maintain the decomposition.
- **Complexity:** Zone theorem gives bounds.

## Idea 9: "Zipper" Algorithm
- **Concept:** If the polygon is a long strip, we can triangulate in one pass.
- **Generalization:** Decompose polygon into "strips" or "serpentines".

## Idea 10: Hybrid Bucket-Splay
- **Concept:** Use buckets to roughly sort/locate, Splay for fine-grained.

## Idea 11: The "Hertel-Mehlhorn" Reverse
- **Concept:** Start with a grid triangulation (which is not valid for the polygon) and fix it? No.

## Selected Approach for Iteration: **Splay-Tree Sweep** & **Bucket-Grid Hybrid**

We will try to implement a **Splay-Tree Sweep**. It is:
1.  **Extremely Simple**: Standard Plane Sweep + Splay Tree.
2.  **Fast**: Likely $O(n)$ for "real" polygons.
3.  **Theoretical Gold Mine**: We can try to prove that the access pattern of a simple polygon on the sweep line satisfies the "Dynamic Finger" property or similar, leading to linear time.

We can also try **Bucket Triangulation** (Uniform Grid) as a robust engineering solution that might be effectively linear.

Let's start with **Splay Sweep**.

