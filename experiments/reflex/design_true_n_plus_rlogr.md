# True O(n + r log r) Algorithm Design

## Why Standard Sweep is O(n log r), Not O(n + r log r)

In standard sweep:
- Every START vertex inserts an edge into sorted status
- Every END vertex deletes an edge
- Every REGULAR_LEFT vertex replaces an edge
- Every SPLIT/MERGE vertex searches for left edge

**Problem**: START/END/REGULAR need to find position in sorted order = O(log r) each
- There are O(n) such vertices
- Total: O(n log r)

## Key Insight for O(n + r log r)

The only operations that NEED O(log r) search are at SPLIT/MERGE vertices.

For all other vertices, we should be able to use O(1) operations via:
1. **Direct pointers**: Know where edges are without searching
2. **Local structure**: Use polygon adjacency for O(1) position finding

## New Algorithm: Monotone-First Decomposition

### Phase 1: Identify Critical Vertices O(n)
Walk boundary, classify all vertices.
Collect r critical vertices (split + merge).

### Phase 2: Sort Critical Vertices O(r log r)  
Sort the r critical vertices by y-coordinate.

### Phase 3: Process Chunks O(n)
The critical vertices divide the sweep into O(r+1) "chunks".
Within each chunk, no split/merge occurs = locally monotone.

**Key**: Within a chunk, we can use simpler stack-based processing without BST!

### Phase 4: Handle Boundaries O(r log r)
At chunk boundaries (critical vertices), we need BST search.
But there are only O(r) such vertices.

## Detailed Algorithm

```
1. Classify all vertices: O(n)
2. Sort critical vertices by y: O(r log r)
3. Initialize status = empty linked list
4. For each chunk from top to bottom:
   a. Process non-critical vertices in chunk using stack: O(chunk size)
   b. At critical vertex, use BST for search: O(log r)
5. Triangulate resulting monotone pieces: O(n)
```

## Why This Works

- Steps 1, 3, 5: O(n)
- Step 2: O(r log r)
- Step 4a: Total O(n) across all chunks (each vertex in exactly one chunk)
- Step 4b: O(r) critical vertices × O(log r) each = O(r log r)

**Total: O(n + r log r)**

## Alternative: Boundary-Order with Reflex-BST

Process vertices in **boundary order** (not y-order):
- Maintain BST of reflex vertices only (size O(r))
- At each vertex, do O(1) local work
- Only reflex vertices interact with BST

This is closer to Hertel-Mehlhorn's approach.

## Most Promising: Lazy Trapezoidalization

1. Build a sparse trapezoidalization using only critical y-levels
2. This requires O(r log r) work
3. Within each trapezoid, triangulate in O(size) time
4. Total: O(n + r log r)

## Implementation Strategy

Let me implement the "Chunk-Based" approach:
1. Sort r critical vertices
2. Process in chunks
3. Use stack for monotone portions, BST only at boundaries

