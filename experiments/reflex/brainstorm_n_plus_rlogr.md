# Brainstorm: True O(n + r log r) Polygon Triangulation

## The Gap

Current algorithm: O(n log r) - all n vertices do O(log r) work
Target: O(n + r log r) - only r vertices do O(log r) work, rest do O(1)

## Key Insight

Which vertices REALLY need O(log r) work?

| Type | Need Search? | Insert/Delete? | Can be O(1)? |
|------|--------------|----------------|--------------|
| START | No | Yes (insert) | Need to find position... |
| END | No | Yes (delete) | YES - have pointer |
| SPLIT | YES | Yes (insert) | NO - must search |
| MERGE | YES | Yes (delete) | NO - must search |
| REGULAR-L | No | Yes (replace) | YES - same position |
| REGULAR-R | Maybe? | No | Depends... |

**Critical observation**: Only SPLIT and MERGE vertices truly need to SEARCH the status to find the left neighbor edge. That's exactly r vertices!

## The Problem with START Vertices

When a START vertex inserts an edge, it needs to find where in the sorted list the edge belongs. This normally requires a search.

But wait... can we avoid this?

## Idea 1: Lazy Landmark Indexing

Maintain two structures:
1. **L**: Doubly-linked list of ALL edges in sorted order
2. **T**: Skip list of O(r) "landmark" edges (sparse index)

**For non-critical vertices (START/END/REGULAR)**:
- Insert/delete in L using local pointer operations
- Don't touch T
- Cost: O(1)

**For critical vertices (SPLIT/MERGE)**:
- Search T to find nearest landmark: O(log r)
- Scan L from that landmark to find exact edge: O(scan)
- Add this edge as new landmark in T: O(log r)
- Total per critical vertex: O(log r + scan)

**Key claim**: Total scan distance across all critical vertices is O(n).

**Proof**: Each edge is scanned at most once. After scanning past an edge, a landmark is placed, preventing future rescans.

**Total time**: O(n) + O(r log r) + O(n) = O(n + r log r)

## Idea 2: Chain-Based Pointers

Observation: The polygon boundary consists of chains. Within a chain, vertices are processed in a specific order.

For a right-chain (vertices where interior is to the left):
- The "left edge" is the same for all vertices in the chain!
- We can track this with a single pointer per chain
- Updates are O(1)

For a left-chain:
- Vertices replace one edge with another
- Position is known from previous edge
- O(1) update

## Idea 3: Precompute Chain Structure

1. **Phase 1**: Walk boundary, identify chains: O(n)
2. **Phase 2**: Sort critical vertices only: O(r log r)
3. **Phase 3**: Process chunks between critical vertices: O(n)
4. **Phase 4**: Handle critical vertices with sparse BST: O(r log r)

## Most Promising: Idea 1 (Lazy Landmarks)

This is the cleanest and most provable approach.

### Algorithm Sketch

```
L = empty doubly-linked list
T = empty skip list
landmarks = {}  // edge -> is_landmark

for v in vertices (sorted by y):
    if v is START:
        e = outgoing edge of v
        insert e into L at position found via local pointers
        
    elif v is END:
        e = incoming edge of v  
        remove e from L
        if e in landmarks: remove from T
        
    elif v is SPLIT:
        // Search for left edge
        landmark = T.find_le(v.x)  // O(log r)
        e_left = scan_right(L, landmark, v.x)  // O(scan)
        if e_left not in landmarks:
            T.insert(e_left)  // O(log r)
            landmarks[e_left] = true
        insert_diagonal(v, helper[e_left])
        helper[e_left] = v
        insert outgoing edge into L
        
    elif v is MERGE:
        // Similar to SPLIT
        ...
        
    elif v is REGULAR:
        // O(1) operations using local pointers
        ...
```

### Why This Works

1. **T has size O(r)**: We add one landmark per split/merge vertex
2. **Searches are O(log r)**: T has O(r) elements
3. **Total scan is O(n)**: Each edge scanned at most once
4. **Non-critical vertices are O(1)**: Only touch L, not T

### Total: O(n + r log r)

## Next Steps

1. Implement this algorithm
2. Verify correctness on test cases
3. Measure operations to confirm O(n + r log r)
4. Write rigorous proofs

