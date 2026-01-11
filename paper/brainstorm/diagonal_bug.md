# Diagonal Generation Bug Analysis

## Observation
For n=200, seed=0 polygon:
- Both chain_optimal and reflex_cli chain_only produce 191 triangles (expected: 198)
- Both generate 57 diagonals
- Both find 56 interior faces (expected: 58 = 57+1)
- Missing 7 triangles = missing 2 faces

## What This Means
The chain-based sweep is generating diagonals that don't properly partition the polygon.
Specifically, some diagonals might:
1. Cross each other (invalid for monotone decomposition)
2. Be placed outside the polygon
3. Connect to wrong helpers

## Root Cause Hypotheses

### H1: Chain advancement timing
When we advance a chain to the sweep line, we might skip over vertices that should
trigger diagonal creation.

### H2: Pending helper management
The pending helper might not be updated correctly when multiple events happen at
similar y-coordinates.

### H3: Left-chain vs right-chain confusion
For "regular on left chain" vertices, we should clear the pending helper, but
we might be checking this incorrectly.

### H4: Numeric precision
Near-collinear points might cause incorrect vertex classification or chain ordering.

## Debugging Plan
1. Visualize the polygon and generated diagonals
2. Check if any diagonals cross
3. Trace the sweep for a specific failing region
4. Compare diagonal generation with a known-correct implementation

## Fallback Options
If the chain sweep is fundamentally flawed:
1. Use edge-based sweep (O(n log n)) for correctness
2. Claim practical speedup from better constants, not asymptotic improvement
3. Document limitations honestly

## Current Status
The chain sweep works for simpler polygons (n <= 100) but fails on some
larger random polygons. This suggests edge cases in the sweep logic rather
than a fundamental algorithmic flaw.
