Dear Editors,

Please consider our manuscript, “Practical polygon triangulation in \(O(n + k\log k)\) time,” for publication in *Computational Geometry: Theory and Applications (CGTA)*.

### Summary and contributions
We present a deterministic triangulation algorithm for simple polygons with \(n\) vertices whose running time is \(O(n + k\log k)\), where \(k\) is the number of local extrema (local maxima, equivalently local minima) with respect to the sweep direction. The algorithm is designed to be both practically fast and implementation-friendly:

- **Instance-sensitive bound in \(k\)**: the sweep processes only \(2k\) extremal events; all regular vertices are handled implicitly by amortized pointer advancement along monotone chains.
- **Chain-based sweep formulation**: we maintain **active monotone chains** rather than individual edges and advance chain pointers lazily, reducing constant factors compared to classical edge-based monotone decomposition.
- **Practical performance**: our implementation is consistently the fastest method in our reported convex, dent, and random benchmarks across all reported sizes, while remaining competitive on highly nonconvex inputs.
- **Reproducibility**: the submission includes a complete reproduction package (`paper/CGAT/`) with scripts to build binaries, run benchmarks, regenerate tables, and rebuild the PDF, along with correctness sanity checks (non-crossing diagonals and \(n-2\) triangles).

### Novelty and relation to prior work
Classical monotone decomposition (Garey et al.) provides an \(O(n\log n)\) triangulation pipeline; Seidel provides an expected-time randomized alternative. Hertel and Mehlhorn also give an instance-sensitive sweep with running time depending on the number of start vertices (an extrema count with respect to the sweep direction). Our work keeps the sweep-based structure but makes the dependence on the **event complexity** \(k\) explicit and implementation-friendly via a **chain-centric status structure**: only the \(2k\) extremal events incur balanced-tree work, while regular vertices are handled in \(O(n)\) total by lazy chain advancement.

### Practical state-of-the-art comparison
We compare against widely used practical implementations: Garey-style monotone decomposition (PolyPartition), Hertel–Mehlhorn via CGAL’s approximate convex partition, and Seidel’s triangulation baseline. The included tables report mean \(\pm\) standard deviation over multiple instances and show clear wall-clock improvements for nearly convex polygons (\(k=O(1)\)) and strong performance on random polygons.

### Availability and reproducibility
All code, benchmark scripts, and generated tables are provided in the repository. The `paper/CGAT/reproduce.sh` script rebuilds the package end-to-end.

Thank you for your time and consideration.

Sincerely,  
Pavel Shpagin, Vasyl Tereschenko
