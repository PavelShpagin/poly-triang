### PDF audit (CGAT submission)

- **Build now succeeds from source**: `paper/CGAT/submission.tex` was missing `algpseudocode` switch/case helpers; added minimal `\\Switch/\\Case/\\EndCase/\\EndSwitch` definitions so `pdflatex` works on TeX Live 2023.
- **Abstract grammar**: fixed “we achieve a about …” → “we achieve about a …”.
- **Convex benchmark block updated**: `paper/CGAT/generated/benchmark_full.tex` convex rows updated to reflect the current optimized code path and no longer show Hertel–Mehlhorn winning on convex.

Build command (WSL):

```bash
cd paper/CGAT
pdflatex -interaction=nonstopmode -halt-on-error submission.tex
pdflatex -interaction=nonstopmode -halt-on-error submission.tex
```

### Code audit (triangulator + harness)

#### Fixed / improved

- **Removed a hidden per-call copy+allocation**: `reflex_tri::Triangulator::triangulate()` (chain version) used to `return triangles_` by value, forcing an extra vector allocation/copy each run. It now returns `const std::vector<Triangle>&` and call sites were updated accordingly.
- **Convex/small‑n fast path**:
  - Removed an unconditional `unordered_set::reserve(...)` that was happening even on convex early-return.
  - Added a small‑n convex pre-check to skip the orientation/reflex passes when the polygon is convex (gated at `n <= 2048` to avoid an extra full scan on large non-convex inputs).
  - Tweaked convex fan construction to use `resize` + indexed writes.
- **Linked implementation**: reserved `triangles` capacity for convex fan path.
- **Benchmark harness ergonomics**: `scripts/benchmark_paper.py` now supports `--types convex,random,...` to avoid running expensive families when you only need a subset.
- **Docs/comments**: updated complexity commentary from `r` to the paper’s `k` parameter where it was misleading.

#### Known follow-ups (not done)

- **Artifact duplication (intentionally kept minimal)**: `paper/CGAT/code/ours/` still duplicates a small subset of the implementation that also exists under `methods/`.
  This is intentional to keep the `paper/CGAT/` artifact self-contained for reproduction. After cleanup, the duplicated surface area is limited to:
  `reflex_chain_triangulate.hpp`, `reflex_fast_linked.hpp`, `reflex_cli.cpp`, and the correctness/debug helper `diag_debug_cli.cpp` (plus a README).
  If we ever decide that self-contained reproduction is not needed, we can fully de-duplicate by building the CGAT binaries directly from `methods/` (but then `paper/CGAT/` would no longer be standalone).
- **Naming collisions**: there are multiple `reflex_tri::Point/Triangle/Triangulator` definitions across different headers; mixing them in one TU will break.
- **Potential perf hot-spot**: face-walk reverse-index construction currently does an \(O(\\sum \\deg(v)^2)\\) scan; if this shows up in profiles, consider building per-vertex “neighbor→slot” maps after angular sort.

