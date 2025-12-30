# True O(n + r log r) Triangulation - Deep Brainstorm

## Why This is Hard

Standard sweep: Every vertex touches BST → O(n log r)
Target: Only r vertices touch BST → O(n) + O(r log r)

## Fundamental Insight

The O(log r) work comes from SEARCHING in a sorted structure.
If we can AVOID searching for non-critical vertices, we win.

---

## Idea 1: Slab Decomposition

**Concept**: Critical y-levels divide plane into r+1 horizontal slabs.
Within each slab, no split/merge occurs → structure is static.

**Algorithm**:
1. Sort r critical y-values: O(r log r)
2. Within each slab, edge order doesn't change
3. Non-critical vertices: O(1) lookup using slab structure
4. Critical vertices: O(log r) BST update

**Problem**: How to assign non-critical vertices to slabs efficiently?

---

## Idea 2: Region Inheritance via Boundary Walk

**Concept**: Walking the boundary, we cross regions in a predictable way.
Each vertex can inherit its region from its boundary predecessor.

**Algorithm**:
1. Precompute region structure at critical levels: O(r log r)
2. Walk boundary CCW, tracking current region: O(n)
3. Each vertex knows its region → O(1) status lookup

**Key Insight**: Boundary walk never "jumps" regions randomly.
We enter/exit regions only at specific edge crossings.

---

## Idea 3: Two-Phase: Skeleton + Fill

**Phase 1** - Reflex Skeleton: O(r log r)
1. Sort reflex vertices by y
2. Process only reflex vertices
3. Build "skeleton" of diagonals among reflex vertices
4. This divides polygon into O(r) "cells"

**Phase 2** - Linear Fill: O(n)
1. Each cell is bounded by skeleton edges + convex chains
2. Triangulate each cell in O(cell size)
3. Sum of cell sizes = n

**Key**: Phase 1 does all the hard work. Phase 2 is trivial.

---

## Idea 4: Event Point Reduction

**Standard sweep**: n events (one per vertex)
**New approach**: r events (critical vertices only)

**Between events**:
- Process vertices in boundary order
- Track position using linked list only (no BST)
- Positions evolve predictably

**At events**:
- Update BST structure: O(log r)
- Adjust linked list pointers

**Problem**: START vertices still need to find position...

---

## Idea 5: START Vertex Special Handling

**Observation**: START vertices don't split regions.
They create new chains but into EXISTING regions.

**Solution**:
- Precompute which region each START vertex falls into
- During sweep, START inserts into known region: O(1)

**How to precompute**:
- At each critical y-level, we know the regions
- For START vertex v between levels y1 and y2:
  - Region at v equals region at y1 or y2 (interpolate)
- Boundary walk tells us which region each START is in

---

## Idea 6: Merge-Sort Vertices Without Full Sort

**Observation**: Vertices on boundary are "almost sorted" by y within chains.

**Algorithm**:
1. Identify r critical vertices: O(n)
2. Sort critical vertices: O(r log r)
3. Non-critical vertices form chains between critical ones
4. Merge chains with sorted critical vertices: O(n)
   - Each chain is monotone, so merge is linear

**Result**: All vertices in y-order in O(n + r log r)

**Then**: Process with position inheritance (Idea 2)

---

## Idea 7: Persistent Chains

**Concept**: Maintain edge information persistently.

1. Process critical vertices, recording "snapshots"
2. Each snapshot has O(r) edges
3. Non-critical vertex v looks up snapshot at v.y: O(log r)
   - But this is still O(log r) per vertex!

**Fix**: Don't lookup. Instead, inherit.

---

## Idea 8: Bucket by Critical Level

1. Create r+1 buckets, one per slab between critical y-levels
2. Each vertex goes into its bucket: O(n)
3. Within each bucket, process in boundary order: O(bucket size)
4. Between buckets (at critical vertices): O(log r)

**Total**: O(n) + O(r log r) = O(n + r log r)

**Key**: Within a bucket, we never need BST because structure is fixed!

---

## Most Promising: Idea 3 (Skeleton + Fill) + Idea 8 (Buckets)

### Algorithm Design

**Phase 0: Classification** O(n)
- Walk boundary, classify all vertices
- Identify r critical vertices

**Phase 1: Bucket Assignment** O(n)
- Sort critical y-values: O(r log r)
- Assign each vertex to a bucket (between two critical levels)
- This is just comparing y-values to sorted critical list

**Phase 2: Process Buckets** O(n + r log r)
- Process buckets top to bottom
- Within bucket: linked list operations only, O(1) per vertex
- At bucket boundary (critical vertex): BST update, O(log r)

### Why It Works

Within a bucket:
- No split/merge vertices exist
- Edge structure is fixed
- Vertices are on monotone chains
- Insert/delete at known positions: O(1)

At bucket boundary:
- Split creates new region: O(log r) to find position
- Merge closes region: O(log r) to find position

### Implementation Details

**Data Structure**:
- Linked list L of edges (sorted by x at current y)
- BST T used ONLY at critical vertices
- Region pointers for each bucket

**Within bucket**:
- START: Insert edge after left boundary of its region
- END: Delete edge (have direct pointer)
- REGULAR: Update in place

**At critical vertex**:
- SPLIT: Search BST for position, update region structure
- MERGE: Search BST for position, merge regions

---

## Let's Implement This!

