# Polygon Triangulation Methods Analysis

## Goal
Achieve O(n + r log r) deterministic polygon triangulation where:
- n = number of vertices
- r = number of monotonicity-violating reflex vertices

## Methods Explored

### 1. Naive Boundary Walking (REJECTED)
- **Approach**: For each reflex vertex, walk the boundary to find support edges
- **Issue**: Same edge can be traversed by multiple reflex vertices
- **Complexity**: O(n * r) = O(n²) worst case
- **Status**: ❌ Not O(n + r log r)

### 2. Persistent Chain Pointers (VERIFIED)
- **Approach**: Process reflex vertices in y-sorted order, splice out regions after each chord
- **Key Insight**: After splicing, traversed edges are removed from the polygon
- **Complexity**: O(n + r log r) via potential function argument
- **Status**: ✅ Correct, implemented and tested
- **Files**: `method_persistent_chains/theory.tex`, `method_persistent_chains/triangulate_v2.py`

### 3. Reflex Graph Method (VERIFIED)
- **Approach**: Build graph G of chord adjacencies, use Reflex-Extrema List for navigation
- **Key Insight**: List has O(r) vertices, shrinks as processing continues
- **Complexity**: O(n + r log r)
- **Status**: ✅ Theoretically sound
- **Files**: `method_reflex_graph/theory.tex`

### 4. Sweep Line without BST (CANCELLED)
- **Approach**: Use sweep line but avoid BST operations
- **Issue**: Hard to avoid O(log n) per query without BST
- **Status**: ❌ Cancelled - doesn't improve on standard approach

## Key Theorems

### Theorem: Non-Intersecting Horizontal Chords
Horizontal chords from distinct reflex vertices do not cross.
- Different heights → parallel lines
- Same height → disjoint segments

### Theorem: O(n) Edge Traversal
Using potential function Φ = edges in remaining polygon:
- Each splice removes ≥ k - O(1) edges (where k = edges traversed)
- Amortized cost per reflex vertex: O(1)
- Total: O(n)

### Theorem: O(n) for Integer Coordinates
With radix sort instead of comparison sort:
- Sort r integers in [0, n^c]: O(r)
- Total: O(n)

## Critical Implementation Details

1. **Linked List Surgery**: After finding chord endpoints, properly disconnect the valley/peak region
2. **Vertex Pointers**: Each reflex vertex needs current pointers to its polygon neighbors
3. **Edge Splitting**: Chord endpoints split polygon edges, creating new vertices

## Recommended Approach

Use **Persistent Chain Pointers** method:
1. Build doubly-linked list from polygon
2. Identify and sort reflex vertices by y-coordinate
3. Process upward vertices bottom-to-top:
   - Walk to find support edges
   - Create chord endpoints
   - Splice out valley
4. Process downward vertices top-to-bottom (symmetric)
5. Triangulate y-monotone regions

## Files Structure

```
experiments/
├── visibility/
│   ├── algorithm.tex     # Main paper (consolidated)
│   └── SUMMARY.md        # This file
├── method_persistent_chains/
│   ├── theory.tex        # Rigorous theory
│   └── triangulate_v2.py # Implementation
└── method_reflex_graph/
    ├── theory.tex        # Alternative theory
    └── triangulate.py    # Implementation (has bugs)
```

## Conclusion

O(n + r log r) is achievable through careful splicing of the polygon structure after each chord addition. The key insight is that edges traversed during support finding are removed from the polygon, ensuring each edge is visited at most once. For integer coordinates, radix sort yields O(n) total time.

