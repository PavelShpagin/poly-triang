# Polygon Triangulation Benchmark Suite

A comprehensive benchmark framework for comparing simple polygon triangulation algorithms.

## Baselines

| CLI | Upstream implementation | Notes |
|-----|------------------------|-------|
| `earcut_cli` | [mapbox/earcut.hpp](https://github.com/mapbox/earcut.hpp) | Canonical ear clipping benchmark |
| `polytri_cli` | [tcoppex/polytri](https://github.com/tcoppex/polytri) | Seidel-style randomized triangulation |
| `cgal_polygon_cli` | [CGAL polygon_triangulation.cpp](https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2polygon_triangulation_8cpp-example.html) | Constrained Delaunay triangulation |

Each CLI simply adapts the upstream example to the `.poly` dataset format and reports runtime.

## Project Structure

```
poly-triang/
├── baselines/          # Algorithm implementations
│   ├── include/        # Headers
│   └── src/            # Source files
├── external/           # Third-party libraries
│   ├── earcut.hpp/     # Mapbox ear clipping
│   └── polytri/        # Seidel's algorithm
├── polygons/           # Test polygon data
│   └── generated/      # Auto-generated test cases
├── results/            # Benchmark output
├── scripts/            # Build and visualization scripts
├── paper/              # LaTeX research paper
└── tests/              # Unit tests
```

## Requirements

- Ubuntu 22.04+ (tested on WSL)
- `sudo apt install build-essential cmake libcgal-dev libcgal-qt5-dev`
- Python 3.10+ with `pandas`, `matplotlib`, `seaborn`

## Quick Start

```bash
./scripts/baseline.sh        # generate polygons, build runners, execute all baselines
python3 scripts/visualize.py # convert CSV → plots + LaTeX table
```

## Output

- `results/benchmark_results.csv` - Raw benchmark data
- `paper/figures/` - Performance plots
- `paper/benchmark_table.tex` - LaTeX table for paper
- `paper/paper.tex` - Complete research paper

## Benchmark matrix

`scripts/generate_polygons.py` produces four deterministic families:

- `convex_n`: regular n-gons (0 reflex vertices)
- `random_n`: random radius star-shaped polygons
- `star_n`: alternating inner/outer radii (≈n/2 reflex vertices)
- `spiral_n`: logarithmic spiral chains

The default sizes are `10, 50, 100, 500, 1000`. Results land in `results/benchmark_results.csv`.
