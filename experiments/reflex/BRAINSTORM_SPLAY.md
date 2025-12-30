# Breakthrough Idea: Self-Adjusting Sweep Line

## The Concept

Instead of a static balanced BST (AVL/Red-Black), use a **Splay Tree** for the sweep-line status.

**Why?**
1. **Dynamic Finger Property**: Splay trees perform well when accesses are spatially close to previous accesses.
2. **Sequential Access Lemma**: Accessing elements in sequential order takes $O(n)$ time.
3. **Working Set Property**: Accessing a small set of items repeatedly is fast.

## Application to Polygons

**Convex Polygon Case ($r=0$):**
- Sweep processes left chain and right chain.
- Left chain edges: $e_{L1}, e_{L2}, \dots$
- Right chain edges: $e_{R1}, e_{R2}, \dots$
- The "active" edges in the BST are just $\{e_L, e_R\}$.
- Updates are strictly sequential/local.
- Splay Tree should handle this in $O(n)$ amortized time.

**General Case:**
- Edges form "bundles" or "chains".
- Within a chain, updates are local.
- Reflex vertices cause "jumps" in the access pattern.
- Splay tree overhead should pay for jumps ($O(\log n)$) but be free for local updates ($O(1)$).

## Expected Complexity

- **Convex**: $O(n)$
- **General**: $O(n + \sum \log \Delta_i)$ where $\Delta_i$ is the "distance" between consecutive operations.
- This effectively captures the $O(n + r \log r)$ or $O(n \log r)$ behavior adaptively.

## Plan

1. Implement a robust Splay Tree.
2. Implement Sweep Line using Splay Tree.
3. Benchmark against Convex and Star polygons.
4. Verify linear behavior on Convex.

