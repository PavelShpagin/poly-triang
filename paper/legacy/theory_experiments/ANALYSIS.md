# Algorithm Variant Analysis

## Why Current Implementation is Slow

### Bottleneck 1: std::set with expensive comparator
- Each BST operation calls `chain_x_at()` for O(log r) pairs
- `chain_x_at()` has lazy advancement loop + interpolation
- Cache misses: each chain accessed in random order

### Bottleneck 2: triangulate_faces() is O(n * deg * log(deg))  
- Line 408-414: For EACH vertex, sort neighbors by angle
- For high-r polygons: many diagonals → high degree → expensive sorting
- This is NOT O(n) as claimed!

### Bottleneck 3: Face extraction with linear search
- Line 453: Linear search for predecessor in neighbor list
- Should use a map or sorted lookup

## Key Insight for Low Constants

For O(n + r log r) to beat O(n log n) Garey when r ≈ n/2:
- Both are O(n log n) asymptotically  
- We need constant factor < 2x Garey's

Garey's strengths:
1. Simple edge comparison (just interpolation)
2. std::set well-optimized
3. No chain overhead

## Variant Ideas

### Variant A: Optimized Edge-Based (Control)
- Standard Garey but with micro-optimizations
- Establishes baseline for what's achievable

### Variant B: Flat Vector with Batched Rebuilds  
- Use sorted vector instead of std::set
- Rebuild every K events instead of incremental updates
- Better cache locality, simpler code

### Variant C: Static Chain Assignment
- Pre-assign each vertex to its chain
- Pre-compute chain x-intercepts at event levels
- O(1) comparisons during sweep

### Variant D: Ear Clipping with Reflex Grid
- O(n²) worst case but O(n) expected
- Extremely simple, cache-friendly
- Skip monotone decomposition entirely

### Variant E: Direct Diagonal Computation
- Compute diagonals without sweep line
- Use polygon structure directly
- Potentially O(n) for many cases

### Variant F: Two-Phase Hybrid
- Phase 1: O(n) identify which vertices need diagonals
- Phase 2: O(r log r) compute diagonal endpoints
- Minimal BST operations

## Success Criteria

1. Beat Garey on convex polygons: 5x+ faster
2. Beat Garey on random CGAL polygons: ANY speedup
3. Never slower than Garey on any polygon type
4. Simple implementation (< 300 lines core algorithm)
5. Rigorous O(n + r log r) proof

