Dear Editors,

Please consider our manuscript, “Practical polygon triangulation in \(O(n + k\log k)\) time,” for publication in *Computational Geometry: Theory and Applications (CGTA)*.

### Summary and contributions
We present a deterministic triangulation algorithm for simple polygons with \(n\) vertices whose running time is \(O(n + k\log k)\), where \(k\) is the number of local extrema (local maxima, equivalently local minima) with respect to the sweep direction. The algorithm is designed to be both practically fast and implementation-friendly:

- **Instance-sensitive bound in \(k\)**: the sweep processes only \(2k\) extremal events; all regular vertices are handled implicitly by amortized pointer advancement along monotone chains.
- **Chain-based sweep formulation**: we maintain **active monotone chains** rather than individual edges and advance chain pointers lazily, reducing constant factors compared to classical edge-based monotone decomposition.
- **Practical performance**: our implementation is consistently the fastest method in our reported convex, dent, and random benchmarks across all reported sizes, while remaining competitive on highly nonconvex inputs.
- **Reproducibility**: the submission includes a self-contained package (`paper/CGAT/`) with scripts to build binaries, run benchmarks, regenerate tables, and rebuild the PDF, along with a correctness sanity check (non-crossing diagonals and \(n-2\) triangles).

### Novelty and relation to prior work
Classical monotone decomposition (Garey et al.) provides an \(O(n\log n)\) triangulation pipeline; Seidel provides an expected-time randomized alternative. Our work keeps the sweep-based structure but makes the complexity and the implementation depend explicitly on the **event complexity** \(k\) rather than on \(n\), and introduces a **chain-centric status structure** that yields an explicit \(O(n)\) amortized handling of regular vertices.

To the best of our knowledge, this explicit extrema-parameterization of monotone decomposition for polygon triangulation—together with an implementation-oriented formulation and empirical evaluation against strong practical baselines—has not been previously stated in this form.

### Practical state-of-the-art comparison
We compare against widely used practical implementations: Garey-style monotone decomposition (PolyPartition), Hertel–Mehlhorn via CGAL’s approximate convex partition, and Seidel’s triangulation baseline. The included tables report mean \(\pm\) standard deviation over multiple instances and show clear wall-clock improvements for nearly convex polygons (\(k=O(1)\)) and strong performance on random polygons.

### Availability and reproducibility
All code, benchmark scripts, and generated tables are provided in the repository. The `paper/CGAT/reproduce.sh` script rebuilds the package end-to-end.

Thank you for your time and consideration.

Sincerely,  
Pavel Shpagin, Vasyl Tereschenko
